#ifndef  bitecoin_endpoint_client_hpp
#define  bitecoin_endpoint_client_hpp

#include <cstdio>
#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include <vector>
#include <memory>
#include <map>
#include <unordered_map>

#include "bitecoin_protocol.hpp"
#include "bitecoin_endpoint.hpp"
#include "bitecoin_hashing.hpp"

#include "cl_utilities.hpp"
#include "cl_miner.hpp"

namespace bitecoin{

struct treeNode
{
    uint32_t index = -1;
    
    treeNode* bitLow = NULL;
    treeNode* bitHigh = NULL;
};
    
class EndpointClient
	: public Endpoint
{
private:
	EndpointClient(EndpointClient &) = delete;
	void operator =(const EndpointClient &) = delete;

	std::string m_minerId, m_clientId;

	unsigned m_knownRounds;
	std::map<std::string,unsigned> m_knownCoins;
    
    //TODO: Worry about zombie memory if you actually go ahead with this:
    
    //        treeNode* root = new treeNode;
    
    //        std::vector<std::vector<uint32_t>> candidates;
    //        std::vector<std::vector<uint32_t>> masterCandidates;
    
    //        for(int i=0; i<160; i++)
    //        {
    //            treeNode* currentNode = root;
    //            for(int j=0;j<256;j++)
    //            {
    //                if (bitIsHigh(rawProofs[i],j))
    //                {
    //                    if (currentNode->bitHigh == NULL) currentNode->bitHigh = new treeNode;
    //                    currentNode = currentNode->bitHigh;
    //                }
    //                else
    //                {
    //                    if (currentNode->bitLow == NULL) currentNode->bitLow = new treeNode;
    //                    currentNode = currentNode->bitLow;
    //                }
    //                if (j==255) currentNode->index = i;
    //            }
    //        }
    
    //findCandidates(candidates, masterCandidates, root);
    
//    void findCandidates(std::vector<std::vector<uint32_t>>& candidateList,std::vector<std::vector<uint32_t>>& masterCandidateList,treeNode* rootNode)
//    {
//        if (rootNode->index != -1)
//        {
//            candidateList.push_back(*new std::vector<uint32_t>(1,rootNode->index));
//        }
//        
//        if (rootNode->bitHigh == NULL && rootNode->bitLow == NULL) return;
//        
//        std::vector<std::vector<uint32_t>> highList;
//        std::vector<std::vector<uint32_t>> lowList;
//        
//        if (rootNode->bitHigh != NULL) findCandidates(highList,masterCandidateList,rootNode->bitHigh);
//        if (rootNode->bitLow != NULL) findCandidates(lowList,masterCandidateList,rootNode->bitLow);
//        
//        if (highList.size() == 0 && lowList.size() == 0) return;
//        
//        std::vector<std::vector<uint32_t>> shortList;
//        
//        shortList.insert(shortList.end(), lowList.begin(), lowList.end());
//        
//        if(highList.size() == 2)
//        {
//            std::vector<uint32_t> combinedHighList(begin(highList[0]),end(highList[0]));
//            combinedHighList.insert(end(combinedHighList), begin(highList[1]),end(highList[1]));
//            shortList.push_back(combinedHighList);
//        }
//        
//        std::sort(std::begin(shortList), std::end(shortList), [] (const std::vector<uint32_t>& left, const std::vector<uint32_t>& right) {
//            return left.size() > right.size();
//        });
//        
//        for(int i=0;i<shortList.size();i++)
//        {
//            if (shortList[0].size() == 16)
//            {
//                masterCandidateList.push_back(shortList[i]);
//            }
//            if (shortList[0].size() < 16)
//            {
//                candidateList.push_back(shortList[i]);
//            }
//        }
//    }
    
    
public:
	
	EndpointClient(
			std::string clientId,
			std::string minerId,
			std::unique_ptr<Connection> &conn,
			std::shared_ptr<ILog> &log
		)
		: Endpoint(conn, log)
		, m_minerId(minerId)
		, m_clientId(clientId)
		, m_knownRounds(0)
	{}
		
	/* Here is a default implementation of make bid.
		I would suggest that you override this method as a starting point.
	*/
	virtual void MakeBid(
		const std::shared_ptr<Packet_ServerBeginRound> roundInfo,	// Information about this particular round
		const std::shared_ptr<Packet_ServerRequestBid> request,		// The specific request we received
		double period,																			// How long this bidding period will last
		double skewEstimate,																// An estimate of the time difference between us and the server (positive -> we are ahead)
		std::vector<uint32_t> &solution,												// Our vector of indices describing the solution
		uint32_t *pProof																		// Will contain the "proof", which is just the value
	){
		double startTime=now()*1e-9;
        
        double tSafetyMargin=0.5;	// accounts for uncertainty in network conditions
		/* This is when the server has said all bids must be produced by, plus the
			adjustment for clock skew, and the safety margin
		*/
		double tFinish=request->timeStampReceiveBids*1e-9 + skewEstimate - tSafetyMargin;
		
		Log(Log_Verbose, "MakeBid - start, total period=%lg.", period);
		
		/*
			We will use this to track the best solution we have created so far.
		*/
		std::vector<uint32_t> bestSolution(roundInfo->maxIndices);
		bigint_t bestProof;
		wide_ones(BIGINT_WORDS, bestProof.limbs);
		
        const Packet_ServerBeginRound *pParams = roundInfo.get();
        
        hash::fnv<64> hasher;
		uint64_t chainHash=hasher((const char*)&pParams->chainData[0], pParams->chainData.size());
        
        //cl_instance instance = create_cl_instance("mining_kernel.cl");
        
        uint32_t microrounds = 100;
        
		unsigned nTrials=0;
		
        std::vector<bigint_t> rawProofs(roundInfo->maxIndices*10);
        
        std::vector<ensemble> candidates;
        
        //TODO: Explore how injecting *greater than the maximum normally reachable random number* can improve score
        
        //TODO: Make sure no constants have been left in here by mistake.
        
        for(unsigned j=0;j<rawProofs.size();j++){
            rawProofs[j] = PoolHash(pParams,j,chainHash);
            
            bigint_t proof = rawProofs[j];
            std::vector<uint32_t> indexes;
            
            indexes.push_back(j);
            
            candidates.push_back(*new ensemble{proof,indexes});
        }
        
        std::sort(std::begin(candidates), std::end(candidates), [] (const ensemble& left, const ensemble& right) {
            return wide_compare(8, left.value.limbs, right.value.limbs) == 1; // Sort in descending order, s.t. smallest value is in highest index
        });
        
        std::vector<ensemble> finalCandidates;
        
        std::vector<uint32_t> usedIndexes;
        
        
        //TODO: Make sure that the zero-hash is never used, as this is unreachable using the random number generator.
        for(int col=255;col>-1;col--)
        {
            int firstNonzeroRow = -1;
            
            for (int row=0;row<160;row++)
            {
                if (bitIsHigh(candidates[row].value,col) && std::find(usedIndexes.begin(), usedIndexes.end(), row) == usedIndexes.end())
                {
                    firstNonzeroRow = row;
                    usedIndexes.push_back(row);
                    break;
                }
            }
            
            if (firstNonzeroRow == -1) continue;
            
            bigint_t firstNonzeroRowValue = candidates[firstNonzeroRow].value;
            std::sort(candidates[firstNonzeroRow].components.begin(),candidates[firstNonzeroRow].components.end());
            
            for (int row=0;row<160;row++)
            {
                if (row == firstNonzeroRow) continue;
                
                if (bitIsHigh(candidates[row].value,col))
                {
                    std::sort(candidates[row].components.begin(),candidates[row].components.end());
                    
                    std::vector<uint32_t> mergedList(candidates[row].components.size()+candidates[firstNonzeroRow].components.size());
                    
                    std::vector<uint32_t>::iterator iterator;
                    
                    iterator=std::set_symmetric_difference(candidates[row].components.begin(),candidates[row].components.end(), candidates[firstNonzeroRow].components.begin(),candidates[firstNonzeroRow].components.end(), mergedList.begin());
                    
                    mergedList.resize(iterator-mergedList.begin());
                    
                    if (mergedList.size() <= 16)
                    {
                        bigint_t tmp;
                        wide_xor(8, tmp.limbs, firstNonzeroRowValue.limbs, candidates[row].value.limbs);
                        candidates[row].value = tmp;
                        candidates[row].components = mergedList;
                    }
                    
                    if (mergedList.size() == 16) finalCandidates.push_back(candidates[row]);
                }
            }
        }
        
        std::sort(std::begin(finalCandidates), std::end(finalCandidates), [] (const ensemble& left, const ensemble& right) {
            return wide_compare(8, left.value.limbs, right.value.limbs) == 1; // Sort in descending order, s.t. smallest value is in highest index
        });
        
        ensemble bestEnsemble = finalCandidates.back();
        
        std::sort(bestEnsemble.components.begin(),bestEnsemble.components.end());
        
        solution=bestEnsemble.components;
        
		wide_copy(BIGINT_WORDS, pProof, bestEnsemble.value.limbs);
        
//        Log(Log_Info, "%lg",wide_as_double(BIGINT_WORDS, rawProofs[0].limbs));
//        Log(Log_Info, "%lg",wide_as_double(BIGINT_WORDS, rawProofs[1].limbs));
//        Log(Log_Info, "%lg",wide_as_double(BIGINT_WORDS, rawProofs[158].limbs));
//        Log(Log_Info, "%lg",wide_as_double(BIGINT_WORDS, rawProofs[159].limbs));
//        
//        Log(Log_Info, "%d",first_zero(rawProofs[159].limbs));
//        Log(Log_Info, "%d",first_zero(rawProofs[158].limbs));
        
//        while(1){
//			++nTrials;
//			
//			Log(Log_Debug, "Trial %d.", nTrials);
//			
//            //bigint_t cl_best_proof = do_gpu_mining(microrounds, pParams, chainHash, instance);
//            
//            std::vector<uint32_t> indices(roundInfo->maxIndices); // Size of maxIndices on test server is 16.
//			uint32_t curr=0;
//			for(unsigned j=0;j<indices.size();j++){
//				curr=curr+1+(rand()%10);
//				indices[j]=curr;
//			}
//			
//			//bigint_t proof=HashReference(pParams, indices.size(), &indices[0],chainHash);
//			
//            bigint_t proof=HashReferenceLite(rawProofs, indices.size(), &indices[0]);
//            double score=wide_as_double(BIGINT_WORDS, proof.limbs);
//			Log(Log_Debug, "    Score=%lg", score);
//			
//			if(wide_compare(BIGINT_WORDS, proof.limbs, bestProof.limbs)<0){
//				double worst=pow(2.0, BIGINT_LENGTH*8);	// This is the worst possible score
//				Log(Log_Verbose, "    Found new best, nTrials=%d, score=%lg, ratio=%lg.", nTrials, score, worst/score);
//				bestSolution=indices;
//				bestProof=proof;
//			}
//			
//			double t=now()*1e-9;	// Work out where we are against the deadline
//			double timeBudget=tFinish-t;
//			Log(Log_Debug, "Finish trial %d, time remaining =%lg seconds.", nTrials, timeBudget);
//			
//			if(timeBudget<=0)
//				break;	// We have run out of time, send what we have
//		}
//		
//		solution=bestSolution;
//		wide_copy(BIGINT_WORDS, pProof, bestProof.limbs);
		
		Log(Log_Verbose, "MakeBid - finish.");
        
        double endTime=now()*1e-9;
        
        Log(Log_Fatal, "Seconds elapsed = %lg seconds.", endTime-startTime);
	}
		
	void Run()
	{
		try{
			auto beginConnect=std::make_shared<Packet_ClientBeginConnect>(m_clientId, m_minerId);
			Log(Log_Info, "Connecting with clientId=%s, minerId=%s", m_clientId.begin(), m_minerId.begin());
			SendPacket(beginConnect);
			
			auto endConnect=RecvPacket<Packet_ServerCompleteConnect>();
			Log(Log_Info, "Connected to exchange=%s, running=%s", endConnect->exchangeId.c_str(), endConnect->serverId.c_str());
			
			while(1){
				Log(Log_Verbose, "Waiting for round to begin.");
				auto beginRound=RecvPacket<Packet_ServerBeginRound>();
				Log(Log_Info, "Round beginning with %u bytes of chain data.", beginRound->chainData.size());
				
				Log(Log_Verbose, "Waiting for request for bid.");
				auto requestBid=RecvPacket<Packet_ServerRequestBid>();
				// Get an estimate of the skew between our clock and theirs. If it is positive,
				// then we are ahead of them.
				double tNow=now()*1e-9;
				double skewEstimate=tNow - requestBid->timeStampRequestBids*1e-9;
				// And work out how long they expect it to last, independent of the skew
				double period=requestBid->timeStampReceiveBids*1e-9 - requestBid->timeStampRequestBids*1e-9;
	
				Log(Log_Info, "Received bid request: serverStart=%lf, ourStart=%lf, skew=%lg. Bid period=%lf", requestBid->timeStampRequestBids*1e-9,  tNow, skewEstimate, period);
				
				std::shared_ptr<Packet_ClientSendBid> bid=std::make_shared<Packet_ClientSendBid>();
				
				MakeBid(beginRound, requestBid, period, skewEstimate, bid->solution, bid->proof);
				bid->timeSent=now();				
				Log(Log_Verbose, "Bid ready.");
				
				SendPacket(bid);
				Log(Log_Verbose, "Bid sent.");
				
				Log(Log_Verbose, "Waiting for results.");
				auto results=RecvPacket<Packet_ServerCompleteRound>();
				Log(Log_Info, "Got round results.");
				
				for(unsigned i=0;i<results->submissions.size();i++){
					double taken=requestBid->timeStampReceiveBids*1e-9 - results->submissions[i].timeRecv*1e-9;
					bool overDue=requestBid->timeStampReceiveBids < results->submissions[i].timeRecv;
					Log(Log_Info, "  %16s : %.6lg, %lg%s", results->submissions[i].clientId.c_str(),
							wide_as_double(BIGINT_WORDS, results->submissions[i].proof), taken,
							overDue?" OVERDUE":""
					);
					if(m_knownCoins.find(results->submissions[i].clientId)==m_knownCoins.end()){
						m_knownCoins[results->submissions[i].clientId]=0;
					}
				}
				
				if(results->winner.clientId==m_clientId){
					Log(Log_Info, "");
					Log(Log_Info, "You won a coin!");
					Log(Log_Info, "");
				}
				
				m_knownRounds++;
				m_knownCoins[results->winner.clientId]++;
				
				Log(Log_Info, "  %16s : %6s, %8s\n", "ClientId", "Coins", "Success");
				auto it=m_knownCoins.begin();
				while(it!=m_knownCoins.end()){
					Log(Log_Info, "  %16s : %6d, %.6lf", it->first.c_str(), it->second, it->second/(double)m_knownRounds);
					++it;
				}
				
				Log(Log_Verbose, "");
			}

		}catch(std::exception &e){
			Log(Log_Fatal, "Exception : %s.", e.what());
			throw;
		}
	}
};

}; // bitecoin

#endif
