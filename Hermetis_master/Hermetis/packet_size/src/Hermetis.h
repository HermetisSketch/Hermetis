#ifndef _LRUQuantile_H
#define _LRUQuantile_H
// #define MEMORY 1

#include "common.h"
#include "leanPart.h"
#include "histogram.h"
#include "coarsePart.h"


template<typename T, int _slot_num, int Interval_num, int _bucket, int _row_number, int _cycle>
class sketch {
public:
    int q_head = 0;
	pair<uint64_t,int> *q = new pair<uint64_t,int> [q_size]{};
	const int count_threshold = COUNT_THRES;
  
	hashTable<T, Interval_num, _slot_num> heavyPart;
	
	Hermetis<_bucket,_row_number,_cycle> lightPart;

	new_data_t flowKeysHist[1000005];
	uint32_t flowKeyPtr = 0;
	uint32_t bandwidth = 0;
	double Time_threshold = (double)_cycle / 3;
	double heavyPart_step = (3.0 * calc_next_prime(_slot_num)) / _cycle;
	double lightPart_step = (3.0 * calc_next_prime(_row_number)) / _cycle;
	
	std::map<k, std::array<T, BUCKET_NUM>> evictHistResult;
	
	std::map<uint64_t, std::array<T, BUCKET_NUM>> cmHistResult;	// decoding result for Histogram

	struct statistics {
		uint32_t pointLarge = 0;
		uint32_t pointAll = 0;
		uint32_t histogramLarge = 0;
		uint32_t histogramAll = 0;
		double pointLargeARE = 0;
		double pointAllARE = 0;
		double histogramLargeARE = 0;
		double histogramAllARE = 0;
		double cardinalityRE = 0;
		double entropyRE = 0;
		uint32_t pointChange = 0;
		uint32_t histogramChange = 0;
	} statistics;

	
    void insert(uint64_t k, uint8_t bid, new_data_t& tmp_key, uint time_tag, unsigned int num, T v = 1) {
		// heavyPart.Clock_Go(num * heavyPart_step);
		// lightPart.Clock_Go(num);

		if((num + 1) % (30000/3) == 0) {
			heavyPart.flush(num);
			lightPart.flush(num);
		}

		uint64_t swap_key;
		slot<T, Interval_num> swap_slot;
		T swap_val = 0;  
        unsigned int swap_num;
		// uint64_t tmpKey(k);
		evictHistResult[tmpKey][bid] += v;
		// cout << " hist check point 0  " << endl;
		int ret = heavyPart.insert(k, bid, num, 30000, swap_key, swap_num, swap_val, swap_slot, v, false);
		// cout << " HIT !!  " << heavyPart.isExist(k) << " " << ((int)lightPart.Query(new_data_t(k, bid).str,num,14)) << " ret == " << ret << endl;
		switch (ret) {
			case 0:{
				uint16_t swap = 0;
				break;
			}
			case HIT_EVICT:				// swap a bucket
			{
				if(((int)num - (int)swap_num ) >= _cycle) break;
				// cout << "swap_val= " << swap_val << endl;
				lightPart.Insert(swap_key, swap_num, swap_val, 14);
				break;
			}
			case MISS_EVICT:			// swap a flow
			{
				for (int i = 0; i < Interval_num; ++i) {
					if ((int)swap_slot.hist.buckets[i].idx == -1) continue;
					new_data_t tmp_key2(swap_slot.key, (uint8_t)swap_slot.hist.buckets[i].idx);
					lightPart.Insert(tmp_key2.str, swap_slot.hist.buckets[i].timestamp[0], swap_slot.hist.buckets[i].Sum(),14);
					bandwidth += 1 + sizeof(T);
				}
				bandwidth += CHARKEY_LEN;
				break;
			}
			case MISS_INSERT:			//just insert into the light part
			{
				lightPart.Insert(tmp_key.str, num, v, 14);
				break;
			}
		}
		return ;	
	}

	int pointQuery(uint64_t k, uint8_t bid, unsigned long long int num, const bool debug = false) {
		uint16_t swap = 0;
		int result = 0;
		// Hash table result
		result = (int)heavyPart.query(k, bid, swap);
		if(heavyPart.isExist(k) == false) result += (int)lightPart.Query(new_data_t(k, bid).str,num,14,debug);
		// if(result == 0) result += (int)lightPart.Query(new_data_t(k, bid).str,num,14,debug);	
		return result;
	}			

	Hist_t histogramQuery(uint64_t key, unsigned long long int num) {
		Hist_t result{}; 
		for (int i = 0; i < BUCKET_NUM; ++i) {
			result.at(i) = pointQuery(key, (uint8_t)i, num);
		}
		return result;
	}

	uint32_t get_cardinality(unsigned long long int num) {
		uint32_t cardinality = lightPart.get_cardinality();
		std::cout << "lightPart.get_cardinality() ==  " << cardinality << std::endl;

		for (int i = 0; i < heavyPart.slot_num; ++i) {
			uint64_t cur_five_tuple = heavyPart.slots[i].key; 
			if (heavyPart.slots[i].hist.getTotal() == 0) {continue;}
				

			Hist_t table_result{};
			std::array<bool, BUCKET_NUM> table_flag{};
			for (int j = 0; j < Interval_num; ++j) {
				if (heavyPart.slots[i].hist.buckets[j].Sum() == 0 )
					continue;
				table_result.at(heavyPart.slots[i].hist.buckets[j].idx) += heavyPart.slots[i].hist.buckets[j].Sum();
			}

			for (int j = 0; j < BUCKET_NUM; ++j) {
				new_data_t cur_key((uint64_t)cur_five_tuple.str, (uint8_t)j);
				T cur_val = table_result.at(j);
				uint32_t light_result = 0;
				// cur_val += lightPart.Query(cur_key.str,num,14);
				light_result = lightPart.Query(cur_key.str,(int)num);

				if (table_flag.at(j) && light_result) {
					cur_val += light_result;
					cardinality--;
				}

				if (cur_val) {
					cardinality++;
				}
			}
		}
		std::cout << "final cardinality ==  " << (int)cardinality << std::endl;
		return cardinality;
	}

	double get_entropy(unsigned long long int num) {
		std::pair<uint32_t, double> entropy = lightPart.get_entropy();
		
		for (int i = 0; i < heavyPart.slot_num; ++i) {
			uint64_t cur_five_tuple = heavyPart.slots[i].key; 
			if (heavyPart.slots[i].hist.getTotal() == 0)
				continue;

			Hist_t table_result{};
			std::array<bool, BUCKET_NUM> table_flag{};

			// initialize
			for (int j = 0; j < BUCKET_NUM; ++j) {
				table_result.at(j) = 0;
			}

			for (int j = 0; j < Interval_num; ++j) {
				if (heavyPart.slots[i].hist.buckets[j].Sum() == 0)
					continue;
				table_result.at(heavyPart.slots[i].hist.buckets[j].idx) += heavyPart.slots[i].hist.buckets[j].Sum();
			}

			for (int j = 0; j < BUCKET_NUM; ++j) {
				new_data_t cur_key((uint64_t)cur_five_tuple.str, (uint8_t)j);
				T cur_val = table_result.at(j);
				uint32_t light_result = lightPart.Query(cur_key.str,(int)num);

				if (table_flag.at(j) && light_result) {
					cur_val += light_result;
					entropy.first -= light_result;
					entropy.second -= light_result * log2(light_result);
				}

				if (cur_val) {
					entropy.first += cur_val;
					entropy.second += cur_val * log2(cur_val);
				}
			}
		}

		return -entropy.second / entropy.first + log2(entropy.first);
	}

	size_t get_memory_usage() {
		std::cout << "\nheavyPart:" << heavyPart.get_memory_usage() << std::endl;
		std::cout << "lightPart: " << lightPart.get_memory_usage() << std::endl;
		return heavyPart.get_memory_usage() + lightPart.get_memory_usage();
	}

};

#endif