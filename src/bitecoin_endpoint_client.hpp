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

#define tbbCores 16
#define tbbOffset 0xFFFFFFFF/tbbCores
#define shortListLengthDefault 3072
#define shortListLengthFast 1024
#define timeGuard 1.7

namespace bitecoin{

    int count_bits_set(uint32_t v)
    {
        // Taken from: http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
        v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
        v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
        return ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count
    }
    
    int wide_hamming_distance(const uint32_t *a, const uint32_t *b)
    {
        bigint_t combined;
        wide_xor(8, combined.limbs, a, b);
        
        uint32_t distance = 0;
        
        for(int i=0;i<8;i++)
        {
            distance+=count_bits_set(combined.limbs[i]);
        }
        return distance;
    }
    
    int wide_hamming_weight(const uint32_t *a)
    {
        uint32_t distance = 0;
        
        distance+=1*count_bits_set(a[4]);
        distance+=5*count_bits_set(a[5]);
        distance+=10*count_bits_set(a[6]);
        distance+=20*count_bits_set(a[7]);
        
        return distance;
    }
    
    int wide_hamming_weight_compare(const uint32_t *a, const uint32_t *b)
    {
        if (a==b) return 0;
        
        uint32_t w_a = wide_hamming_weight(a);
        uint32_t w_b = wide_hamming_weight(b);
        
        if (w_a < w_b) return -1;
        if (w_a > w_b) return 1;
        return 0;
    }
    
	class EndpointClient
		: public Endpoint
	{
	private:
		EndpointClient(EndpointClient &) = delete;
		void operator =(const EndpointClient &) = delete;

		std::string m_minerId, m_clientId;

		unsigned m_knownRounds;
		std::map<std::string, unsigned> m_knownCoins;

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
			double startTime = now()*1e-9;
			double tSafetyMargin = timeGuard;
			double tFinish = request->timeStampReceiveBids*1e-9 + skewEstimate - tSafetyMargin;

			Log(Log_Verbose, "MakeBid - start, total period=%lg.", period);

			const Packet_ServerBeginRound *pParams = roundInfo.get();

			// This doesn't change within each round, so calculate it once and re-use.
			hash::fnv<64> hasher;
			uint64_t chainHash = hasher((const char*)&pParams->chainData[0], pParams->chainData.size());

			std::vector<ensemble> candidates;
            
            double t = now()*1e-9;
            double timeBudget = tFinish - t;
            
            uint32_t shortListLength = timeBudget > 1.0 ? shortListLengthDefault : shortListLengthFast;
            
			auto compMin = [](const ensemble& left, const ensemble& right) {
				return wide_compare(8,left.value.limbs, right.value.limbs) == 1;
			};

			auto compMax = [](const ensemble& left, const ensemble& right) {
				return wide_compare(8,left.value.limbs, right.value.limbs) == -1;
			};

			tbb::task_group group;

			std::vector<std::priority_queue<ensemble, std::vector<ensemble>, decltype(compMin)>> priorityQueues;
			std::vector<uint32_t> totalTrials(tbbCores);

			for (int i = 0; i < tbbCores; i++)
			{
				std::priority_queue<ensemble, std::vector<ensemble>, decltype(compMin)> ensemble_priority_queue(compMin);

				priorityQueues.push_back(ensemble_priority_queue);

				group.run([&, i](){

					std::priority_queue<ensemble, std::vector<ensemble>, decltype(compMax)> ensemble_priority_queue_reversed(compMax);

					unsigned nTrials = 0;
					unsigned offset = tbbOffset * i;
					while (1)
					{
						bigint_t proof = PoolHash(pParams, nTrials + offset, chainHash);

						if (priorityQueues[i].size() < shortListLength || wide_compare(8,proof.limbs, ensemble_priority_queue_reversed.top().value.limbs) == -1)
						{
							std::vector<uint32_t> indexes;

							indexes.push_back(nTrials + offset);

							ensemble e = ensemble{ proof, indexes };

							priorityQueues[i].push(e);
							ensemble_priority_queue_reversed.push(e);
						}

						if (ensemble_priority_queue_reversed.size() > shortListLength) ensemble_priority_queue_reversed.pop();

						double t = now()*1e-9;
						double timeBudget = tFinish - t;

						Log(Log_Debug, "Finish trial %d, time remaining =%lg seconds.", nTrials, timeBudget);

						nTrials++;

						if ((timeBudget <= 0 && priorityQueues[i].size() >= shortListLength) || nTrials >= tbbOffset-1)
						{
							totalTrials[i] = nTrials;
							break;	// We have run out of time, send what we have
						}

					}
				});

			}

			group.wait();

			uint32_t overallTrials = std::accumulate(totalTrials.begin(), totalTrials.end(), 0);

			for (int i = 0; i < shortListLength; i++)
			{
				auto nextQueue = std::min_element(priorityQueues.begin(), priorityQueues.end(), [](const std::priority_queue<ensemble, std::vector<ensemble>, decltype(compMin)>& left, const std::priority_queue<ensemble, std::vector<ensemble>, decltype(compMin)>& right) { return wide_compare(8, left.top().value.limbs, right.top().value.limbs) == -1;
				});

				candidates.push_back(nextQueue->top());
				nextQueue->pop();
			}

			Log(Log_Info, "Tried %d elements", overallTrials);

			double gStart = now()*1e-9;
            
            for(int r=0;r<4;r++)
            {
                uint32_t c_size = candidates.size();
                
                std::vector<ensemble> newCandidates;
                
                tbb::task_group group;
                
                std::vector<std::priority_queue<ensemble, std::vector<ensemble>, decltype(compMin)>> priorityQueues;
                
                for (int t = 0; t < tbbCores; t++)
                {
                    std::priority_queue<ensemble, std::vector<ensemble>, decltype(compMin)> ensemble_priority_queue(compMin);
                    
                    priorityQueues.push_back(ensemble_priority_queue);
                    
                    group.run([&, t](){
                        
                        std::priority_queue<ensemble, std::vector<ensemble>, decltype(compMax)> ensemble_priority_queue_reversed(compMax);
                        
                        for(uint32_t i=0;i<c_size/tbbCores;i++)
                        {
                            uint32_t rowA = t*(c_size/tbbCores) + i;
                            
                            for(uint32_t rowB=rowA+1;rowB<c_size;rowB++)
                            {
                                std::vector<uint32_t> mergedList(candidates[rowA].components.size() + candidates[rowB].components.size());
                                
                                std::vector<uint32_t>::iterator iterator;
                                
                                iterator = std::set_symmetric_difference(candidates[rowA].components.begin(), candidates[rowA].components.end(), candidates[rowB].components.begin(), candidates[rowB].components.end(), mergedList.begin());
                                
                                mergedList.resize(iterator - mergedList.begin());
                                
                                if (mergedList.size()==0) continue;
                                
                                ensemble e;
                                wide_xor(8,e.value.limbs,candidates[rowA].value.limbs,candidates[rowB].value.limbs);
                                
                                e.components = mergedList;
                                
                                if (priorityQueues[t].size() < c_size || wide_compare(8,e.value.limbs, ensemble_priority_queue_reversed.top().value.limbs) == -1)
                                {
                                    priorityQueues[t].push(e);
                                    ensemble_priority_queue_reversed.push(e);
                                }
                                
                                if (ensemble_priority_queue_reversed.size() > c_size) ensemble_priority_queue_reversed.pop();

                            }
                        }
                        
                    });
                }
                
                group.wait();
                
                for (int i = 0; i < c_size; i++)
                {
                    auto nextQueue = std::min_element(priorityQueues.begin(), priorityQueues.end(), [](const std::priority_queue<ensemble, std::vector<ensemble>, decltype(compMin)>& left, const std::priority_queue<ensemble, std::vector<ensemble>, decltype(compMin)>& right) { return wide_compare(8, left.top().value.limbs, right.top().value.limbs) == -1;
                    });
                    
                    newCandidates.push_back(nextQueue->top());
                    nextQueue->pop();
                }
                
                candidates = newCandidates;
                
            }

			// Sort the finalists in descending order, s.t. smallest value is in highest index
			std::sort(std::begin(candidates), std::end(candidates), [](const ensemble& left, const ensemble& right) {
				return wide_compare(8, left.value.limbs, right.value.limbs) == 1;
			});

            ensemble bestEnsemble = candidates.back();

            std::sort(bestEnsemble.components.begin(), bestEnsemble.components.end());

            solution = bestEnsemble.components;

            wide_copy(BIGINT_WORDS, pProof, bestEnsemble.value.limbs);

            double endTime = now()*1e-9;

            Log(Log_Verbose, "MakeBid - finish.");
			Log(Log_Info, "Pair Testing Time Used = %lg seconds.", endTime - gStart);
			Log(Log_Info, "Overall Time Used = %lg seconds.", endTime - startTime);
		}

		void Run()
		{
			try{
				auto beginConnect = std::make_shared<Packet_ClientBeginConnect>(m_clientId, m_minerId);
				Log(Log_Info, "Connecting with clientId=%s, minerId=%s", m_clientId.begin(), m_minerId.begin());
				SendPacket(beginConnect);

				auto endConnect = RecvPacket<Packet_ServerCompleteConnect>();
				Log(Log_Info, "Connected to exchange=%s, running=%s", endConnect->exchangeId.c_str(), endConnect->serverId.c_str());

				while (1){
					Log(Log_Verbose, "Waiting for round to begin.");
					auto beginRound = RecvPacket<Packet_ServerBeginRound>();
					Log(Log_Info, "Round beginning with %u bytes of chain data.", beginRound->chainData.size());

					Log(Log_Verbose, "Waiting for request for bid.");
					auto requestBid = RecvPacket<Packet_ServerRequestBid>();
					// Get an estimate of the skew between our clock and theirs. If it is positive,
					// then we are ahead of them.
					double tNow = now()*1e-9;
					double skewEstimate = tNow - requestBid->timeStampRequestBids*1e-9;
					// And work out how long they expect it to last, independent of the skew
					double period = requestBid->timeStampReceiveBids*1e-9 - requestBid->timeStampRequestBids*1e-9;

					Log(Log_Info, "Received bid request: serverStart=%lf, ourStart=%lf, skew=%lg. Bid period=%lf", requestBid->timeStampRequestBids*1e-9, tNow, skewEstimate, period);

					std::shared_ptr<Packet_ClientSendBid> bid = std::make_shared<Packet_ClientSendBid>();

					// Making the bid
					MakeBid(beginRound, requestBid, period, skewEstimate, bid->solution, bid->proof);
					bid->timeSent = now();
					Log(Log_Verbose, "Bid ready.");

					SendPacket(bid);
					Log(Log_Verbose, "Bid sent.");

					Log(Log_Verbose, "Waiting for results.");
					auto results = RecvPacket<Packet_ServerCompleteRound>();
					Log(Log_Info, "Got round results.");

					for (unsigned i = 0; i < results->submissions.size(); i++){
						double taken = requestBid->timeStampReceiveBids*1e-9 - results->submissions[i].timeRecv*1e-9;
						bool overDue = requestBid->timeStampReceiveBids < results->submissions[i].timeRecv;
						Log(Log_Info, "  %16s : %.6lg, %lg%s", results->submissions[i].clientId.c_str(),
							wide_as_double(BIGINT_WORDS, results->submissions[i].proof), taken,
							overDue ? " OVERDUE" : ""
							);
						if (m_knownCoins.find(results->submissions[i].clientId) == m_knownCoins.end()){
							m_knownCoins[results->submissions[i].clientId] = 0;
						}
					}

					if (results->winner.clientId == m_clientId){
						Log(Log_Info, "");
						Log(Log_Info, "You won a coin!");
						Log(Log_Info, "");
					}

					m_knownRounds++;
					m_knownCoins[results->winner.clientId]++;

					Log(Log_Info, "  %16s : %6s, %8s\n", "ClientId", "Coins", "Success");
					auto it = m_knownCoins.begin();
					while (it != m_knownCoins.end()){
						Log(Log_Info, "  %16s : %6d, %.6lf", it->first.c_str(), it->second, it->second / (double)m_knownRounds);
						++it;
					}

					Log(Log_Verbose, "");
				}

			}
			catch (std::exception &e){
				Log(Log_Fatal, "Exception : %s.", e.what());
				throw;
			}
		}
	};

}; // bitecoin

#endif