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

//#include "tbb/task_group.h"
#include "tbb/tbb.h"

#define tbbCores 8
#define tbbOffset 0xFFFFFFFF/tbbCores
#define shortListLengthDefault 2000
#define shortListLengthFast 2000
#define timeGuard 1.5

namespace bitecoin{

	int count_bits_set(uint32_t v)
	{
		v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
		v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
		return ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count
	}

	int wide_hamming_distance(const uint32_t *a, const uint32_t *b)
	{
		bigint_t combined;
		wide_xor(8, combined.limbs, a, b);

		uint32_t distance = 0;

		for (int i = 0; i < 8; i++)
		{
			distance += count_bits_set(combined.limbs[i]);
		}
		return distance;
	}

	int wide_hamming_weight(const uint32_t *a)
	{
		uint32_t distance = 0;

		distance += 1 * count_bits_set(a[4]);
		distance += 5 * count_bits_set(a[5]);
		distance += 10 * count_bits_set(a[6]);
		distance += 20 * count_bits_set(a[7]);

		return distance;
	}

	int wide_hamming_weight_compare(const uint32_t *a, const uint32_t *b)
	{
		if (a == b) return 0;

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
				return wide_compare(8, left.value.limbs, right.value.limbs) == 1;
			};

			auto compMax = [](const ensemble& left, const ensemble& right) {
				return wide_compare(8, left.value.limbs, right.value.limbs) == -1;
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

						if (priorityQueues[i].size() < shortListLength || wide_compare(8, proof.limbs, ensemble_priority_queue_reversed.top().value.limbs) == -1)
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

						if ((timeBudget <= 0 && priorityQueues[i].size() >= shortListLength) || nTrials >= tbbOffset - 1)
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

			//std::reverse(candidates.begin(),candidates.end());

			// This is where we store all the best combinations of xor'ed vectors. Each combination is of size roundInfo->maxIndices
			std::vector<ensemble> finalCandidates;

			for (int r = 0; r < 4; r++)
			{
				std::vector<std::vector<ensemble>> newCandidates(tbbCores, std::vector<ensemble>(candidates.capacity() / tbbCores));

				// Originally for(int i=0;i<candidates.size();i++)
				tbb::parallel_for(0, tbbCores, 1, [&](unsigned x)
				{
					//Reclaim original vector pointer
					std::vector<ensemble>* const newCandidates_ptr = &newCandidates[x];

					for (int z = 0; z < candidates.size()/tbbCores; z++)
					{
						//Reclaim original iterator
						const unsigned i = z + x*candidates.size()/tbbCores;

						for (int j = i + 1; j < candidates.size(); j++)
						{
							ensemble e;
							wide_xor(8, e.value.limbs, candidates[i].value.limbs, candidates[j].value.limbs);

							std::vector<uint32_t> mergedList(candidates[i].components.size() + candidates[j].components.size());

							std::vector<uint32_t>::iterator iterator;

							iterator = std::set_symmetric_difference(candidates[i].components.begin(), candidates[i].components.end(), candidates[j].components.begin(), candidates[j].components.end(), mergedList.begin());

							mergedList.resize(iterator - mergedList.begin());

							if (mergedList.size() == 0) continue;

							e.components = mergedList;

							(*newCandidates_ptr).push_back(e);
						}
					}
				});

				// Now merge newCandidates
				//std::vector<ensemble> newCandidates_merge(candidates.size);
				std::vector<ensemble> newCandidates_merge(candidates.capacity());
				for (unsigned k = 0; k < tbbCores; k++)
					newCandidates_merge.insert(newCandidates_merge.end(), newCandidates[k].begin(), newCandidates[k].end());

				// Sort vector
				std::sort(std::begin(newCandidates_merge), std::end(newCandidates_merge), [](const ensemble& left, const ensemble& right) {
					return wide_hamming_weight_compare(left.value.limbs, right.value.limbs) == -1;
				});

				// Output best
				newCandidates_merge.resize(candidates.size());
				candidates = newCandidates_merge;
			}



			// We find optimal combinations of the proofs calculated for each index using 'Gaussian elimination' (but xor-ing instead of adding/subtracting). We start in the column of the MSB, and xor vectors that have this bit high together to make the bit in this column 0 for as many vectors as possible. We then move to the next most significant bit, and repeat the process. At each stage, we keep track of what set of indexes we are xor-ing with what orther set of indexes. The the combined set size reaches roundInfo->maxIndices, we add this candidate solution to finalCandidates.

			//			//// -- Gaussian Elimination Starts Here -- ////
			//
			//			std::vector<uint32_t> usedIndexes;
			//
			//			for (int col = 255; col > -1; col--)
			//			{
			//				int firstNonzeroRow = -1;
			//
			//				for (int row = 0; row < candidates.size(); row++)
			//				{
			//					if (bitIsHigh(candidates[row].value, col) && std::find(usedIndexes.begin(), usedIndexes.end(), row) == usedIndexes.end())
			//					{
			//						firstNonzeroRow = row;
			//						usedIndexes.push_back(row);
			//						break;
			//					}
			//				}
			//
			//				if (firstNonzeroRow == -1) continue;
			//
			//				bigint_t firstNonzeroRowValue = candidates[firstNonzeroRow].value;
			//				std::sort(candidates[firstNonzeroRow].components.begin(), candidates[firstNonzeroRow].components.end());
			//
			//				for (int row = 0; row < candidates.size(); row++)
			//				{
			//					if (row == firstNonzeroRow) continue;
			//
			//					if (bitIsHigh(candidates[row].value, col))
			//					{
			//						std::sort(candidates[row].components.begin(), candidates[row].components.end());
			//
			//						std::vector<uint32_t> mergedList(candidates[row].components.size() + candidates[firstNonzeroRow].components.size());
			//
			//						std::vector<uint32_t>::iterator iterator;
			//
			//						iterator = std::set_symmetric_difference(candidates[row].components.begin(), candidates[row].components.end(), candidates[firstNonzeroRow].components.begin(), candidates[firstNonzeroRow].components.end(), mergedList.begin());
			//
			//						mergedList.resize(iterator - mergedList.begin());
			//
			//						if (mergedList.size() <= roundInfo->maxIndices)
			//						{
			//							bigint_t tmp;
			//							wide_xor(8, tmp.limbs, firstNonzeroRowValue.limbs, candidates[row].value.limbs);
			//							candidates[row].value = tmp;
			//							candidates[row].components = mergedList;
			//						}
			//
			//						if (mergedList.size() == roundInfo->maxIndices) finalCandidates.push_back(candidates[row]);
			//					}
			//				}
			//			}
			//
			//			//// -- Gaussian Elimination Ends Here -- ////

			// Sort the finalists in descending order, s.t. smallest value is in highest index
			std::sort(std::begin(candidates), std::end(candidates), [](const ensemble& left, const ensemble& right) {
				return wide_compare(8, left.value.limbs, right.value.limbs) == 1;
			});

			// Choose the finalist with the lowest score for our final bid.
			if (!candidates.empty())
			{
				ensemble bestEnsemble = candidates.back();

				std::sort(bestEnsemble.components.begin(), bestEnsemble.components.end());

				solution = bestEnsemble.components;

				wide_copy(BIGINT_WORDS, pProof, bestEnsemble.value.limbs);
			}
			else
			{
				// Last ditch attempt to make sure we always submit something valid. Ideally we should never come in here.

				std::vector<uint32_t> indices(roundInfo->maxIndices);
				uint32_t curr = 0;
				for (unsigned j = 0; j < indices.size(); j++){
					curr = curr + 1 + (rand() % 10);
					indices[j] = curr;
				}

				bigint_t proof = HashReference(pParams, indices.size(), &indices[0], chainHash);

				solution = indices;

				wide_copy(BIGINT_WORDS, pProof, proof.limbs);
			}

			double gEnd = now()*1e-9;

			Log(Log_Info, "GE Time Elapsed = %lg seconds.", gEnd - gStart);

			Log(Log_Verbose, "MakeBid - finish.");

			double endTime = now()*1e-9;

			Log(Log_Info, "Time used = %lg seconds.", endTime - startTime);
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