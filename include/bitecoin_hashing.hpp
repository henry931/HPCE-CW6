#ifndef bitecoin_hashing_hpp
#define bitecoin_hashing_hpp

// Provies a basic (non-cryptographic) hash function
#include "contrib/fnv.hpp"

// Basic operations for dealing with multi-word integers.
// The individual words are often called "limbs", so a 256 bit
// integer contains 8 x 32-bit limbs.
#include "wide_int.h"

#include <cstdint>
#include <stdexcept>
#include <algorithm>
#include <random>

#include "bitecoin_protocol.hpp"

namespace bitecoin{
	
	// The size of the multi-word types we use in hashing
	enum{ NLIMBS = BIGINT_WORDS };
	
	// Convenience wrapper around a multi-word type
	struct bigint_t
	{
		uint32_t limbs[NLIMBS];
	};
	
	// This provides a primitive randomness step. It is not cryptographic quality,
	// but suffices for these purposes. There is a constant c that comes from the
	// server at the beginning of the round that gets used here.
	void PoolHashStep(bigint_t &x, const Packet_ServerBeginRound *pParams)
	{
		assert(NLIMBS==4*2);
		
		bigint_t tmp;
		// tmp=lo(x)*c;
		wide_mul(4, tmp.limbs+4, tmp.limbs, x.limbs, pParams->c);
		// [carry,lo(x)] = lo(tmp)+hi(x)
		uint32_t carry=wide_add(4, x.limbs, tmp.limbs, x.limbs+4);
		// hi(x) = hi(tmp) + carry
            wide_add(4, x.limbs+4, tmp.limbs+4, carry);
		
		// overall:  tmp=lo(x)*c; x=tmp>hi(x)
	}
	
	// Given the various round parameters, this calculates the hash for a particular index value.
	// Multiple hashes of different indices will be combined to produce the overall result.
	bigint_t PoolHash(const Packet_ServerBeginRound *pParams, uint32_t index,uint64_t chainHash)
	{
		assert(NLIMBS==4*2);
        
		// Incorporate the existing block chain data - in a real system this is the
		// list of transactions we are signing. This is the FNV hash:
		// http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
		
		// The value x is 8 words long (8*32 bits in total)
		// We build (MSB to LSB) as  [ chainHash ; roundSalt ; roundId ; index ]
		
        //TODO: Can the 6 lines below be done with one SSE/AVX packing command?
        
        bigint_t x;
		wide_zero(8, x.limbs);
		wide_add(8, x.limbs, x.limbs, index);	//chosen index goes in at two low limbs
		wide_add(6, x.limbs+2, x.limbs+2, pParams->roundId);	// Round goes in at limbs 3 and 2
		wide_add(4, x.limbs+4, x.limbs+4, pParams->roundSalt);	// Salt goes in at limbs 5 and 4
		wide_add(2, x.limbs+6, x.limbs+6, chainHash);	// chainHash at limbs 7 and 6
		
		// Now step forward by the number specified by the server
		for(unsigned j=0;j<pParams->hashSteps;j++){
			PoolHashStep(x, pParams);
		}
		return x;
	}	
	
    struct ensemble
    {
        bigint_t value;
        std::vector<uint32_t> components;
    };
    
    bool bitIsHigh(bigint_t& number, uint32_t position)
    {
        uint32_t word = position/32;
        return number.limbs[word] & (0x1 << (position % 32));
    }
    
    uint32_t nlz1(uint32_t x)
    {
        int n;
        
        if (x == 0) return(-1);
        n = 0;
        if (x <= 0x0000FFFF) {n = n +16; x = x <<16;}
        if (x <= 0x00FFFFFF) {n = n + 8; x = x << 8;}
        if (x <= 0x0FFFFFFF) {n = n + 4; x = x << 4;}
        if (x <= 0x3FFFFFFF) {n = n + 2; x = x << 2;}
        if (x <= 0x7FFFFFFF) {n = n + 1;}
        return 32-n;
    }
    
    uint32_t first_zero(const uint32_t *limbs)
    {
        for(int i=0;i<8;i++)
        {
            if (limbs[7-i]==0) continue;
            int position = nlz1(limbs[7-i]);
            return position + (7-i)*32;
        }
        return -1;
    }
    
    
	bigint_t HashReferenceLite(std::vector<bigint_t>& rawProofs, unsigned nIndices, const uint32_t *pIndices)
    {
        bigint_t acc;
		wide_zero(8, acc.limbs);
		
		for(unsigned i=0;i<nIndices;i++){
			
			// Combine the hashes of the points together using xor
			wide_xor(8, acc.limbs, acc.limbs, rawProofs[pIndices[i]].limbs);
		}
		
		return acc;

    }
    
    
    // This is the complete hash reference function. Given the current round parameters,
	// and the solution, which is a vector of indices, it calculates the proof. The goodness
	// of the solution is determined by the numerical value of the proof.
	bigint_t HashReference(
		const Packet_ServerBeginRound *pParams,
		unsigned nIndices,
		const uint32_t *pIndices,
        const uint64_t chainHash
	){
		if(nIndices>pParams->maxIndices)
			throw std::invalid_argument("HashReference - Too many indices for parameter set.");
		
		bigint_t acc;
		wide_zero(8, acc.limbs);
		
		for(unsigned i=0;i<nIndices;i++){
			if(i>0){
				if(pIndices[i-1] >= pIndices[i])
					throw std::invalid_argument("HashReference - Indices are not in monotonically increasing order.");
			}
			
			// Calculate the hash for this specific point
			bigint_t point=PoolHash(pParams, pIndices[i],chainHash);
			
			// Combine the hashes of the points together using xor
			wide_xor(8, acc.limbs, acc.limbs, point.limbs);
		}
		
		return acc;
	}
	
	/*! This is used to choose the winner. It is somewhat biased against the very fastest
		people, so that it is still possible for slow people to occasionally win a coin.
		\param rng returns an double-precision uniform random in [0,1)
		\note It is basically choosing in proportion to the sqrt of the inverse fitness.
			So if there are three people, with effective hash rates of 1000, 100, and 10,
			then the probability of each winning a coin is 7%, 22%, and 70%
	*/
	template<class TRng>
	unsigned ChooseWinner(const std::vector<bigint_t> &choices, TRng &rng)
	{
		std::vector<std::pair<double,unsigned> > values;
		for(unsigned i=0;i<choices.size();i++){
			values.push_back(std::make_pair(1.0 / wide_as_double(8, choices[i].limbs), i));
		}
		std::sort(values.begin(), values.end());
		double acc=0;
		for(unsigned i=0;i<values.size();i++){
			acc+=sqrt(values[i].first);
			values[i].first=acc;
		}
		
		std::uniform_real_distribution<double> distribution(0.0,acc);
		
		double sel=distribution(rng);
		for(unsigned i=0;i<values.size();i++){
			if(sel<=values[i].first)
				return values[i].second;
		}
		return values.back().second;
	}
	
};

#endif
