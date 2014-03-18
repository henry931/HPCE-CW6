#ifndef wide_int_h
#define wide_int_h

#include <stdint.h>
#include <assert.h>

#include <stdio.h>
#include <math.h>

// AVX
#include "immintrin.h"

/*! Simply library for maintaining large positive integers as an array
	of 32-bit limbs */

/*! Compare two integers as numbers:
		a<b :  -1
		a==b :  0
		a>b : +1
*/
int wide_compare(unsigned n, const uint32_t *a, const uint32_t *b)
{
	if(a==b)
		return 0;

	for(int i=n-1;i>=0;i--){
		if(a[i]<b[i])
			return -1;
		if(a[i]>b[i])
			return +1;
	}
	return 0;
}

/*! Copy a source number to a destination */
void wide_copy(unsigned n, uint32_t *res, const uint32_t *a)
{
	for(unsigned i=0;i<n;i++){
		res[i]=a[i];
	}
}

/*! Set entire number to zero */
void wide_zero(unsigned n, uint32_t *res)
{
	for(unsigned i=0;i<n;i++){
		res[i]=0;
	}
}

/*! Set entire number to zero */
void wide_ones(unsigned n, uint32_t *res)
{
	for(unsigned i=0;i<n;i++){
		res[i]=0xFFFFFFFFul;
	}
}

/*! Add together two n-limb numbers, returning the carry limb.
	\note the output can also be one of the inputs
*/
void wide_xor(unsigned n, uint32_t *res, const uint32_t *a, const uint32_t *b)
{
	for(unsigned i=0;i<n;i++){
		res[i]=a[i]^b[i];
	}
}

/*! Add together two n-limb numbers, returning the carry limb.
	\note the output can also be one of the inputs
*/
uint32_t wide_add(unsigned n, uint32_t *res, const uint32_t *a, const uint32_t *b)
{
	uint64_t carry=0;
	for(unsigned i=0;i<n;i++){
		uint64_t tmp=uint64_t(a[i])+b[i]+carry;
		res[i]=uint32_t(tmp&0xFFFFFFFFULL);
		carry=tmp>>32;
	}
	return carry;
}


	/*! Add a single limb to an n-limb number, returning the carry limb
	\note the output can also be the input
*/
uint32_t wide_add(unsigned n, uint32_t *res, const uint32_t *a, uint32_t b)
{
	uint64_t carry=b;
	for(unsigned i=0;i<n;i++){
		uint64_t tmp=a[i]+carry;
		res[i]=uint32_t(tmp&0xFFFFFFFFULL);
		carry=tmp>>32;
	}
	return carry;
}

uint32_t wide_add(unsigned n, uint32_t *res, const uint32_t *a, uint64_t b)
{
	assert(n>=2);
	uint64_t acc=a[0]+(b&0xFFFFFFFFULL);
	res[0]=uint32_t(acc&0xFFFFFFFFULL);
	uint64_t carry=acc>>32;
	
	acc=a[1]+(b&0xFFFFFFFFULL)+carry;
	res[1]=uint32_t(acc&0xFFFFFFFFULL);
	carry=acc>>32;
	
	for(unsigned i=2;i<n;i++){
		uint64_t tmp=a[i]+carry;
		res[i]=uint32_t(tmp&0xFFFFFFFFULL);
		carry=tmp>>32;
	}
	return carry;
}

/*! Multiply two n-limb numbers to produce a 2n-limb result
	\note All the integers must be distinct, the output cannot overlap the input */
void wide_mul(unsigned n, uint32_t *res_hi, uint32_t *res_lo, const uint32_t *a, const uint32_t *b)
{
	assert(res_hi!=a && res_hi!=b);
	assert(res_lo!=a && res_lo!=b);
	
	uint64_t carry=0, acc=0;
	for(unsigned i=0; i<n; i++){
		for(unsigned j=0; j<=i; j++){
			assert( (j+(i-j))==i );
			uint64_t tmp=uint64_t(a[j])*b[i-j];
			acc+=tmp;
			if(acc < tmp)
				carry++;
			//fprintf(stderr, " (%d,%d)", j,i-j);
		}
		res_lo[i]=uint32_t(acc&0xFFFFFFFFull);
		//fprintf(stderr, "\n  %d : %u\n", i, res_lo[i]);
		acc= (carry<<32) | (acc>>32);
		carry=carry>>32;
	}
	
	for(unsigned i=1; i<n; i++){
		for(unsigned j=i; j<n; j++){
			uint64_t tmp=uint64_t(a[j])*b[n-j+i-1];
			acc+=tmp;
			if(acc < tmp)
				carry++;
			//fprintf(stderr, " (%d,%d)", j,n-j+i-1);
			//assert( (j+(n-j))==n+i );
		}
		res_hi[i-1]=uint32_t(acc&0xFFFFFFFFull);
		//fprintf(stderr, "\n  %d : %u\n", i+n-1, res_hi[i-1]);
		acc= (carry<<32) | (acc>>32);
		carry=carry>>32;
	}
	res_hi[n-1]=acc;
}

//! Return x as a double, which is obviously lossy for large n
double wide_as_double(unsigned n, const uint32_t *x)
{
	double acc=0;
	for(unsigned i=0;i<n;i++){
		acc+=ldexp((double)x[i], i*32);
	}
	return acc;
}

void wide_mul_avx(uint32_t *res_hi, uint32_t *res_lo, const uint32_t *a, const uint32_t *b)
{
	// So the idea is to do long multiplication...

	// Pack 32-bit limbs
	const __m128i ainput = _mm_set_epi32(a[3], a[2], a[1], a[0]);
	const __m128i binput = _mm_set_epi32(b[3], b[2], b[1], b[0]);
	const __m128i aupper = _mm_cvtepu32_epi64(_mm_srli_si128(ainput, 8));
	const __m128i alower = _mm_cvtepu32_epi64(ainput);
	const __m128i bupper_shift0 = _mm_cvtepu32_epi64(_mm_srli_si128(binput, 8));	// also blower_shift2
	const __m128i blower_shift0 = _mm_cvtepu32_epi64(binput);						// also bupper_shift2
	const __m128i bupper_shift1 = _mm_cvtepu32_epi64(_mm_shuffle_epi32(binput, _MM_SHUFFLE(0, 0, 0, 3))); // also blower_shift3
	const __m128i blower_shift1 = _mm_cvtepu32_epi64(_mm_shuffle_epi32(binput, _MM_SHUFFLE(0, 0, 2, 1))); // also bupper_shift3

	// Multiply limbs
	const __m128i abupper_shift0 = _mm_mul_epu32(aupper, bupper_shift0);
	const __m128i ablower_shift0 = _mm_mul_epu32(alower, blower_shift0);
	const __m128i abupper_shift1 = _mm_mul_epu32(aupper, bupper_shift1);
	const __m128i ablower_shift1 = _mm_mul_epu32(alower, blower_shift1);
	const __m128i abupper_shift2 = _mm_mul_epu32(aupper, blower_shift0);
	const __m128i ablower_shift2 = _mm_mul_epu32(alower, bupper_shift0);
	const __m128i abupper_shift3 = _mm_mul_epu32(aupper, blower_shift1);
	const __m128i ablower_shift3 = _mm_mul_epu32(alower, bupper_shift1);

	// Convert to 64-bit for addition overflow
	const __m128i r0c0 = _mm_cvtepu32_epi64(ablower_shift0);
	const __m128i r0c1 = _mm_cvtepu32_epi64(ablower_shift1);
	const __m128i r0c2 = _mm_cvtepu32_epi64(ablower_shift2);
	const __m128i r2c2 = _mm_cvtepu32_epi64(abupper_shift2);
	const __m128i r0c3 = _mm_cvtepu32_epi64(ablower_shift3);
	const __m128i r2c3 = _mm_cvtepu32_epi64(abupper_shift3);
	const __m128i r2c4 = _mm_cvtepu32_epi64(abupper_shift0);
	const __m128i r2c5 = _mm_cvtepu32_epi64(abupper_shift1);

	// Need to shift these guys 8-bytes right first for the rest
	const __m128i r1c1 = _mm_cvtepu32_epi64(_mm_srli_si128(ablower_shift3, 8));
	const __m128i r1c2 = _mm_cvtepu32_epi64(_mm_srli_si128(ablower_shift0, 8));
	const __m128i r1c3 = _mm_cvtepu32_epi64(_mm_srli_si128(ablower_shift1, 8));
	const __m128i r3c3 = _mm_cvtepu32_epi64(_mm_srli_si128(abupper_shift1, 8));
	const __m128i r1c4 = _mm_cvtepu32_epi64(_mm_srli_si128(ablower_shift2, 8));
	const __m128i r3c4 = _mm_cvtepu32_epi64(_mm_srli_si128(abupper_shift2, 8));
	const __m128i r3c5 = _mm_cvtepu32_epi64(_mm_srli_si128(abupper_shift3, 8));
	const __m128i r3c6 = _mm_cvtepu32_epi64(_mm_srli_si128(abupper_shift0, 8));

	// OK now add
	__m128i result0 = r0c0;
	__m128i result1 = _mm_add_epi64(r0c1, r1c1);
	__m128i result2 = _mm_add_epi64(_mm_add_epi64(r0c2, r1c2), r2c2);
	__m128i result3 = _mm_add_epi64(_mm_add_epi64(r0c3, r1c3), _mm_add_epi64(r2c3, r3c3));
	__m128i result4 = _mm_add_epi64(_mm_add_epi64(r1c4, r2c4), r3c4);
	__m128i result5 = _mm_add_epi64(r2c5, r3c5);
	__m128i result6 = r3c6;

	// Pointers for access
	uint32_t* const res0 = (uint32_t*)&result0;
	uint32_t* const res1 = (uint32_t*)&result1;
	uint32_t* const res2 = (uint32_t*)&result2;
	uint32_t* const res3 = (uint32_t*)&result3;
	uint32_t* const res4 = (uint32_t*)&result4;
	uint32_t* const res5 = (uint32_t*)&result5;
	uint32_t* const res6 = (uint32_t*)&result6;

	uint64_t* const res0_64 = (uint64_t*)&result0;
	uint64_t* const res1_64 = (uint64_t*)&result1;
	uint64_t* const res2_64 = (uint64_t*)&result2;
	uint64_t* const res3_64 = (uint64_t*)&result3;
	uint64_t* const res4_64 = (uint64_t*)&result4;
	uint64_t* const res5_64 = (uint64_t*)&result5;
	uint64_t* const res6_64 = (uint64_t*)&result6;

	// Push overflows into subsequent buffers
	res1_64[0] += res0_64[1] + res0[1];
	res2_64[0] += res1_64[1] + res1[1];
	res3_64[0] += res2_64[1] + res2[1];
	res4_64[0] += res3_64[1] + res3[1];
	res5_64[0] += res4_64[1] + res4[1];
	res6_64[0] += res5_64[1] + res5[1];

	// Output lowest 32-bit in each buffer
	res_lo[0] = res0[0];
	res_lo[1] = res1[0];
	res_lo[2] = res2[0];
	res_lo[3] = res3[0];
	res_hi[0] = res4[0];
	res_hi[1] = res5[0];
	res_hi[2] = res6[0];
	res_hi[3] = res6[2] + res6[1];

}

#endif
