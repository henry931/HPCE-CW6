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
#include <queue>
#include <numeric>

#include "bitecoin_protocol.hpp"
#include "bitecoin_endpoint.hpp"
#include "bitecoin_hashing.hpp"

#include "tbb/task_group.h"

namespace bitecoin{
    
    class EndpointClient
	: public Endpoint
    {
    private:
        EndpointClient(EndpointClient &) = delete;
        void operator =(const EndpointClient &) = delete;
        
        std::string m_minerId, m_clientId;
        
        unsigned m_knownRounds;
        std::map<std::string,unsigned> m_knownCoins;
        
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
        
        virtual void MakeBid(
                             const std::shared_ptr<Packet_ServerBeginRound> roundInfo,
                             const std::shared_ptr<Packet_ServerRequestBid> request,
                             double period,
                             double skewEstimate,
                             std::vector<uint32_t> &solution,
                             uint32_t *pProof
                             ){
            double startTime=now()*1e-9;
            double tSafetyMargin=1.0;
            double tFinish=request->timeStampReceiveBids*1e-9 + skewEstimate - tSafetyMargin;
            
            Log(Log_Verbose, "MakeBid - start, total period=%lg.", period);
            
            const Packet_ServerBeginRound *pParams = roundInfo.get();
            
            // This doesn't change within each round, so calculate it once and re-use.
            hash::fnv<64> hasher;
            uint64_t chainHash=hasher((const char*)&pParams->chainData[0], pParams->chainData.size());
            
            std::vector<ensemble> candidates;
            
            //TODO: Make sure no constants have been left in here by mistake.
            
            auto compMin = [] (const ensemble& left, const ensemble& right) {
                return wide_compare(8, left.value.limbs, right.value.limbs) == 1;
            };
            
            auto compMax = [] (const ensemble& left, const ensemble& right) {
                return wide_compare(8, left.value.limbs, right.value.limbs) == -1;
            };
            
            tbb::task_group group;
            
            std::vector<std::priority_queue<ensemble,std::vector<ensemble>,decltype(compMin)>> priorityQueues;
            std::vector<uint32_t> totalTrials(8);
            
            for(int i=0;i<8;i++)
            {
                std::priority_queue<ensemble,std::vector<ensemble>,decltype(compMin)> ensemble_priority_queue(compMin);
                
                priorityQueues.push_back(ensemble_priority_queue);
                
                group.run([&, i](){
                    
                    std::priority_queue<ensemble,std::vector<ensemble>,decltype(compMax)> ensemble_priority_queue_reversed(compMax);
                    
                    unsigned nTrials = 0;
                    unsigned offset = 2000000*i;
                    while(1)
                    {
                        bigint_t proof = PoolHash(pParams,nTrials+offset,chainHash);
                        
                        if (priorityQueues[i].size()<2048 || wide_compare(8, proof.limbs, ensemble_priority_queue_reversed.top().value.limbs) == -1)
                        {
                            std::vector<uint32_t> indexes;
                            
                            indexes.push_back(nTrials+offset);
                            
                            ensemble* e = new ensemble{proof,indexes};
                            
                            priorityQueues[i].push(*e);
                            ensemble_priority_queue_reversed.push(*e);
                        }
                        
                        if (ensemble_priority_queue_reversed.size() > 2048) ensemble_priority_queue_reversed.pop();
                        
                        double t=now()*1e-9;
                        double timeBudget=tFinish-t;
                        
                        Log(Log_Debug, "Finish trial %d, time remaining =%lg seconds.", nTrials, timeBudget);
                        
                        nTrials++;
                        
                        if(timeBudget<=0 && priorityQueues[i].size() >= 2048)
                        {
                            totalTrials[i] = nTrials;
                            break;	// We have run out of time, send what we have
                        }
                        
                    }
                });
                
            }
            
			group.wait();
            
            uint32_t overallTrials = std::accumulate(totalTrials.begin(),totalTrials.end(),0);
  
            for (int i=0; i<2000; i++)
            {
                auto nextQueue = std::min_element(priorityQueues.begin(),priorityQueues.end(),[] (const std::priority_queue<ensemble,std::vector<ensemble>,decltype(compMin)>& left, const std::priority_queue<ensemble,std::vector<ensemble>,decltype(compMin)>& right) { return wide_compare(8, left.top().value.limbs, right.top().value.limbs) == -1;
                });
                
                candidates.push_back(nextQueue->top());
                nextQueue->pop();
            }
            
            Log(Log_Info, "Tried %d elements", overallTrials);
            
            double gStart=now()*1e-9;
            
            // This is where we store all the best combinations of xor'ed vectors. Each combination is of size roundInfo->maxIndices
            std::vector<ensemble> finalCandidates;
            
            // We find optimal combinations of the proofs calculated for each index using 'Gaussian elimination' (but xor-ing instead of adding/subtracting). We start in the column of the MSB, and xor vectors that have this bit high together to make the bit in this column 0 for as many vectors as possible. We then move to the next most significant bit, and repeat the process. At each stage, we keep track of what set of indexes we are xor-ing with what orther set of indexes. The the combined set size reaches roundInfo->maxIndices, we add this candidate solution to finalCandidates.
            
            //// -- Gaussian Elimination Starts Here -- ////
            
            std::vector<uint32_t> usedIndexes;
            
            for(int col=255;col>-1;col--)
            {
                int firstNonzeroRow = -1;
                
                for (int row=0;row<candidates.size();row++)
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
                
                for (int row=0;row<candidates.size();row++)
                {
                    if (row == firstNonzeroRow) continue;
                    
                    if (bitIsHigh(candidates[row].value,col))
                    {
                        std::sort(candidates[row].components.begin(),candidates[row].components.end());
                        
                        std::vector<uint32_t> mergedList(candidates[row].components.size()+candidates[firstNonzeroRow].components.size());
                        
                        std::vector<uint32_t>::iterator iterator;
                        
                        iterator=std::set_symmetric_difference(candidates[row].components.begin(),candidates[row].components.end(), candidates[firstNonzeroRow].components.begin(),candidates[firstNonzeroRow].components.end(), mergedList.begin());
                        
                        mergedList.resize(iterator-mergedList.begin());
                        
                        if (mergedList.size() <= roundInfo->maxIndices)
                        {
                            bigint_t tmp;
                            wide_xor(8, tmp.limbs, firstNonzeroRowValue.limbs, candidates[row].value.limbs);
                            candidates[row].value = tmp;
                            candidates[row].components = mergedList;
                        }
                        
                        if (mergedList.size() == roundInfo->maxIndices) finalCandidates.push_back(candidates[row]);
                    }
                }
            }
            
            //// -- Gaussian Elimination Ends Here -- ////
            
            // Sort the finalists in descending order, s.t. smallest value is in highest index
            std::sort(std::begin(finalCandidates), std::end(finalCandidates), [] (const ensemble& left, const ensemble& right) {
                return wide_compare(8, left.value.limbs, right.value.limbs) == 1;
            });
            
            // Choose the finalist with the lowest score for our final bid.
            if (!finalCandidates.empty())
            {
                ensemble bestEnsemble = finalCandidates.back();
                
                std::sort(bestEnsemble.components.begin(),bestEnsemble.components.end());
                
                solution=bestEnsemble.components;
                
                wide_copy(BIGINT_WORDS, pProof, bestEnsemble.value.limbs);
            }
            else
            {
                // Last ditch attempt to make sure we always submit something valid. Ideally we should never come in here.
                
                std::vector<uint32_t> indices(roundInfo->maxIndices);
                uint32_t curr=0;
                for(unsigned j=0;j<indices.size();j++){
                    curr=curr+1+(rand()%10);
                    indices[j]=curr;
                }
                
                bigint_t proof=HashReference(pParams, indices.size(), &indices[0],chainHash);
                
                solution = indices;
                
                wide_copy(BIGINT_WORDS, pProof, proof.limbs);
            }
            
            double gEnd=now()*1e-9;
            
            Log(Log_Info, "GE Time Elapsed = %lg seconds.", gEnd-gStart);
            
            Log(Log_Verbose, "MakeBid - finish.");
            
            double endTime=now()*1e-9;
            
            Log(Log_Info, "Time used = %lg seconds.", endTime-startTime);
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