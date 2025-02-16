#ifndef _HEAVY_H
#define _HEAVY_H

#include "histogram.h"
#include "common.h"
#include "utils.h"

#define HT_b (1.08)

template<typename T, int bucket_num>
class slot {
public:
	uint64_t key;					// flowkey
	histogram<T, bucket_num> hist;	// histogram
	T negativeCounter;				// negative counter
	unsigned int T_first,T_last;

	// slot(uint64_t k, uint8_t bid, T v = 0) {
	// 	key = k;
	// 	hist.insert(bid, v);
	// 	negativeCounter = 0;
	// }

	slot(uint64_t k) {
		key = k;
		negativeCounter = 0;
		T_first = 0;
		T_last = 0;
	}

	slot() {
		negativeCounter = 0;
		T_first = 0;
		T_last = 0;
	}

	static size_t get_memory_usage() {
		return sizeof(T) + CHARKEY_LEN + histogram<T, bucket_num>::get_memory_usage();
	}
};

template <typename T, int bucket_num, int size>
class hashTable {
public:
	/* key + histogram + negative_vote */
	slot<T, bucket_num> *slots;
	int nonempty;
	uint32_t slot_num;
	uint64_t h, s, n;
	unsigned int last_time,clock_pos;
	int cycle_num,field_num;

	hashTable() {
		int index = 0;
		h = GenHashSeed(index++);
		s = GenHashSeed(index++);
		n = GenHashSeed(index++);
		nonempty = 0;
		slot_num = calc_next_prime(size);
		slots = new slot<T, bucket_num>[slot_num];
		last_time = 0;
		clock_pos = 0;
		cycle_num = 0;
		field_num = 4;
	}

	uint32_t calPos(uint64_t key) {
		return (uint32_t)(AwareHash((unsigned char *)key, CHARKEY_LEN, h, s, n) % slot_num);
	}

	int insert(uint64_t k, uint8_t bid, unsigned int num, int cycle, new_data_t &swap_key, unsigned int &swap_num, T &swap_val, slot<T, bucket_num> &swap_slot,
				uint32_t v = 1, bool debug = false)					// v: number of packets / pakcet length
	{
		uint64_t tempKey = k;
		uint32_t pos = calPos(k);
		bucket<T> b;

		cycle_num = ((num) / (30000 / 3)) % 4;
		// if((num + 1) % (30000/3) == 0) flush(num);

		if (slots[pos].key == tempKey) {
			int ret = slots[pos].hist.insert(bid, num, cycle, v, b, (cycle_num ) % field_num, false);
			if (ret == HIST_EVICT) {									// evict a bucket of the histogram
				swap_key = new_data_t(k, b.idx);
				
				swap_num = b.timestamp[0];
				swap_val = 0;
				for(int j=0; j<3; ++j){
					swap_val += b.val[j];
				}
				
				return HIT_EVICT;
			}
			
			return HIT;
		}
		else if (slots[pos].key.isEmpty() || slots[pos].hist.getTotal() == 0) {
			slots[pos].key = tempKey;
			slots[pos].hist.insert(bid, num, cycle, v, b, (cycle_num ) % field_num, false);
			nonempty++;
			// if (debug)
			// 	std::cout << "HIT (an empty slot) " << pos << std::endl;
			return HIT;
		}
				T total = slots[pos].hist.getTotal();
				slots[pos].negativeCounter += v;
				if (!(rand()%int(pow(HT_b,(double)(total - slots[pos].negativeCounter))))){
					swap_slot = slots[pos];
					unsigned long long int temp = slots[pos].negativeCounter;
					slots[pos].hist.reset();
					slots[pos] = slot<T, bucket_num>(k);
					slots[pos].negativeCounter = 0;
					slots[pos].hist.insert(bid, num, cycle, v, b, (cycle_num ) % field_num);
					return MISS_EVICT;
				}
				
			
		}

		// if (debug)
		// 	std::cout << "MISS_INSERT\n";

		return MISS_INSERT;
	}

	T query(uint64_t k, uint8_t bid, uint16_t &swap) {
		uint32_t pos = calPos(k);
		swap = 0;

		if (tempKey == slots[pos].key) {
			T result = slots[pos].hist.query(bid, swap);
			return result;
		}
		else {
			return 0;
		}
	}

	bool isExist(uint64_t k){
		uint64_t tempKey = k;
		uint32_t pos = calPos(k);
		if (tempKey == slots[pos].key) {
			return true;
		}
		return false;
	}

	void Clock_Go(unsigned int num){
		for(;last_time < num; ++last_time){
			for(int j = 0;j<bucket_num;++j){
				slots[clock_pos].hist.buckets[j].val[(cycle_num + 1) % field_num] = 0;
			}
			clock_pos = (clock_pos + 1) % slot_num;
			if(clock_pos == 0){
				cycle_num = (cycle_num + 1) % field_num;
			}
		}
	}


	void flush(unsigned int num){
		int interval = ((num) / (30000 / 3) ) % 4;
		int _cycle = 30000;
		for(int i = 0;i < slot_num;++i){
			for(int j = 0;j<bucket_num;++j){
				slots[i].hist.buckets[j].val[(interval + 1) % field_num] = 0;
				
			}
		}
	}

	double calLoadRate() {
		return (double)nonempty / slot_num;
	}

	size_t get_memory_usage() {
		return slot_num * slot<T, bucket_num>::get_memory_usage();
	}
};

#endif