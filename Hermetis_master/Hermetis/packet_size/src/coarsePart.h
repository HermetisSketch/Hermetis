#ifndef _Hermetis_H
#define _Hermetis_H

#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <cstring>
#include <string.h>
#include <stdlib.h>
#include "BOBHash.h"
#include "AwareHash.h"
#include <stdint.h>
#define NDEBUG
#include <cassert>
#include <unordered_map>
#include<bitset>
#include "hash_class.h"
#include <unordered_set>

#define ENABLE_CM_SKETCH 0
/* the bucket type macro */
#define FINGERPRINT_LENGTH 8
#define TIMETAG_LENGTH 2
#define TYPE_ID_LENGTH 3
#define MAXLOOP 3

typedef struct {
	uint8_t type_id:4;
	uint8_t skip:4;
	uint8_t skip2[2];
	uint8_t fingerprint5;
	uint8_t fingerprint4;
	uint8_t fingerprint3;
	uint8_t fingerprint2;
	uint8_t fingerprint1;
} ec_bucket_type_origin;
#define define_type_5_fingerprint(type_index, c1, c2, c3, c4, c5)  typedef struct { \
		uint64_t type_id: TYPE_ID_LENGTH; \
		uint64_t count1: c1; \
		uint64_t count2: c2; \
		uint64_t count3: c3; \
		uint64_t count4: c4; \
		uint64_t count5: c5; \
		uint64_t fingerprint1: FINGERPRINT_LENGTH; \
		uint64_t fingerprint2: FINGERPRINT_LENGTH; \
		uint64_t fingerprint3: FINGERPRINT_LENGTH; \
		uint64_t fingerprint4: FINGERPRINT_LENGTH; \
		uint64_t fingerprint5: FINGERPRINT_LENGTH; \
	} ec_bucket_type##type_index;

#define define_type_4_fingerprint(type_index, c1, c2, c3, c4)  typedef struct { \
		uint64_t type_id: TYPE_ID_LENGTH; \
		uint64_t count1: c1; \
		uint64_t count2: c2; \
		uint64_t count3: c3; \
		uint64_t count4: c4; \
		uint64_t fingerprint1: FINGERPRINT_LENGTH; \
		uint64_t fingerprint2: FINGERPRINT_LENGTH; \
		uint64_t fingerprint3: FINGERPRINT_LENGTH; \
		uint64_t fingerprint4: FINGERPRINT_LENGTH; \
	} ec_bucket_type##type_index;

#define define_type_3_fingerprint(type_index, c1, c2, c3)  typedef struct { \
		uint64_t type_id: TYPE_ID_LENGTH; \
		uint64_t count1: c1; \
		uint64_t count2: c2; \
		uint64_t count3: c3; \
		uint64_t fingerprint1: FINGERPRINT_LENGTH; \
		uint64_t fingerprint2: FINGERPRINT_LENGTH; \
		uint64_t fingerprint3: FINGERPRINT_LENGTH; \
	} ec_bucket_type##type_index;

/* Each type struct */
// Store the bucket item size, each BUCKET_ITEM_SIZE[i] is a type of bucket
// BUCKET_ITEM_SIZE[0]:
//       0: bucket_type_index, i.e. the type_id of the bucket
//       5: bucket_fingerprint_size, i.e. how many fingerprints are in the bucket
//       2, 3, 4, 5, 6: the each item bit size
#define BUCKET_TYPE_NUM 7
#define FINGERPRINT_MAX_NUM 4
#define BUCKET_LENGTH 1
#define CONFIG_LENGTH (2+FINGERPRINT_MAX_NUM)
const uint8_t BUCKET_ITEM_SIZE[BUCKET_TYPE_NUM][CONFIG_LENGTH] = {
	{0, 4, 3, 5, 6, 7},
	{1, 3, 5, 6, 20, 0},
	{2, 3, 6, 7, 18, 0},
	{3, 3, 7, 8, 16, 0},
	{4, 3, 8, 9, 14, 0},
	{5, 3, 9, 10, 12, 0},
	{6, 3, 10, 10, 11, 0}
};

define_type_5_fingerprint(0, 2, 3, 4, 5, 6);
typedef union {
	uint64_t value;
	ec_bucket_type0 type_0;
	ec_bucket_type_origin type0;
} ec_bucket;

struct Unit{
	ec_bucket* count;
};

unordered_set<uint> hotKey_Pos_st;

/* The item localation and  */
// const uint64_t TYPE_ID_LOC = 0; 
// const uint64_t TYPE_ID_BITMASK = (1UL << TYPE_ID_LENGTH) - 1;
#define TYPE_ID_LOC (0)
#define TYPE_ID_BITMASK ((1UL << TYPE_ID_LENGTH) - 1)
// fingerprint_loc: the start bit index of the fingerprint, [x, x+FINGERPRINT_LENGTH]
// fingerprint_loc_mask: the mask for the fingerprint
int FINGERPRINT_LOC[FINGERPRINT_MAX_NUM];
// uint64_t FINGERPRINT_BITMASK = (1UL << FINGERPRINT_LENGTH) - 1; // The mask is 0xFF (8 bits)
#define FINGERPRINT_BITMASK ((1UL << FINGERPRINT_LENGTH) - 1)
#define TIMETAG_BITMASK ((1UL << TIMETAG_LENGTH) - 1)
#define BUCKET_BITMASK ((1UL << BUCKET_LENGTH) - 1)

uint8_t COUNT_LEN[BUCKET_TYPE_NUM][FINGERPRINT_MAX_NUM];
uint8_t _align1[4];
uint8_t TIMETAG_LOC[BUCKET_TYPE_NUM][FINGERPRINT_MAX_NUM];
uint8_t COUNT_LOC[BUCKET_TYPE_NUM][FINGERPRINT_MAX_NUM];
uint8_t _align2[4];

#define get_bucket_type(bucket) (bucket.type0.type_id)
#define set_bucket_type(bucket, type) ((bucket).type0.type_id = (type))

static inline void init_bucket_parameters() {
	// Check the config is correct!
	assert(sizeof(ec_bucket) == 8);
	for (int i = 0; i < BUCKET_TYPE_NUM; i++) {
		const int fingerprint_num = BUCKET_ITEM_SIZE[i][1];
		int total_bit = 0;
		total_bit += TYPE_ID_LENGTH;
		total_bit += FINGERPRINT_LENGTH * fingerprint_num;
		for (int j = 2; j < 2+fingerprint_num; j++) {
			total_bit += 2 + BUCKET_ITEM_SIZE[i][j];
		}
		assert(total_bit == 8 * sizeof(ec_bucket));
		for (int j = 2+fingerprint_num; j < CONFIG_LENGTH; j++) {
			assert(BUCKET_ITEM_SIZE[i][j] == 0);
		}
	}

	for(int i = 0; i < FINGERPRINT_MAX_NUM; i++) {
		FINGERPRINT_LOC[i] = i * FINGERPRINT_LENGTH + TYPE_ID_LENGTH;
		// FINGERPRINT_LOC[i] = 64 - (i+1) * FINGERPRINT_LENGTH; // 0-4, from the highest to smallest
	}

	for(int i = 0; i < BUCKET_TYPE_NUM; i++) {
		const int fingerprint_num = BUCKET_ITEM_SIZE[i][1];
		for (int j = 2; j < 2 + fingerprint_num; j++) {
			const int count_idx = j - 2;
			COUNT_LEN[i][count_idx] = BUCKET_ITEM_SIZE[i][j];
			if ( count_idx == 0 ) {
				TIMETAG_LOC[i][count_idx] = TYPE_ID_LENGTH + FINGERPRINT_LENGTH * fingerprint_num;
				COUNT_LOC[i][count_idx] = TYPE_ID_LENGTH + FINGERPRINT_LENGTH * fingerprint_num + TIMETAG_LENGTH * fingerprint_num;
			} else {
				TIMETAG_LOC[i][count_idx] = TIMETAG_LOC[i][count_idx-1] + TIMETAG_LENGTH;
				COUNT_LOC[i][count_idx] = COUNT_LOC[i][count_idx-1] + COUNT_LEN[i][count_idx-1];
			}
		}
	}
}

static inline __attribute__((always_inline)) uint32_t get_bucket_type_id(ec_bucket* bkt) {
	return (bkt->value & (TYPE_ID_BITMASK));
}

static inline __attribute__((always_inline)) void set_bucket_type_id(ec_bucket* bkt, uint64_t type_id) {
	bkt->value = (bkt->value & ~(TYPE_ID_BITMASK)) | (type_id & TYPE_ID_BITMASK);
}

static inline __attribute__((always_inline)) uint8_t get_item_num_in_bucket_type(uint32_t type_id) {
	// return BUCKET_ITEM_SIZE[type_id][1];
	switch(type_id) {
		case 0: return 4;
		case 1 ... 6: return 3;
		default: return 3;
	}
}

static inline __attribute__((always_inline)) uint8_t get_bucket_fingerprint(ec_bucket* bkt, int fingerprint_index) {
	int temp = (bkt->value >> FINGERPRINT_LOC[fingerprint_index]);
	uint8_t result = ((bkt->value >> FINGERPRINT_LOC[fingerprint_index]) & FINGERPRINT_BITMASK);
	return ((bkt->value >> FINGERPRINT_LOC[fingerprint_index]) & FINGERPRINT_BITMASK);
	// switch (fingerprint_index) {
	// 	case 0: return (bkt->type0.fingerprint1);
	// 	case 1: return (bkt->type0.fingerprint2);
	// 	case 2: return (bkt->type0.fingerprint3);
	// 	case 3: return (bkt->type0.fingerprint4);
	// 	case 4: return (bkt->type0.fingerprint5);
	// 	default: assert(0);
	// }
}

static inline __attribute__((always_inline)) void set_bucket_fingerprint(ec_bucket* bkt, int fingerprint_index, uint64_t fingerprint) {
	bkt->value = (bkt->value & ~(FINGERPRINT_BITMASK << FINGERPRINT_LOC[fingerprint_index])) | ((fingerprint & FINGERPRINT_BITMASK) << FINGERPRINT_LOC[fingerprint_index]);
}

static inline __attribute__((always_inline)) uint64_t get_bucket_time_tag(ec_bucket* bkt, int time_tag_index, const uint32_t type_id) {
	return (bkt->value >> TIMETAG_LOC[type_id][time_tag_index]) & TIMETAG_BITMASK;
}

static inline __attribute__((always_inline)) void set_bucket_time_tag(ec_bucket* bkt, int time_tag_index, uint64_t tag_value, const uint32_t type_id) {
	bkt->value = (bkt->value & ~((TIMETAG_BITMASK) << TIMETAG_LOC[type_id][time_tag_index])) | ((tag_value & (TIMETAG_BITMASK)) << TIMETAG_LOC[type_id][time_tag_index]);
}

static inline __attribute__((always_inline)) uint64_t get_bucket_count(ec_bucket* bkt, int count_index, const uint32_t type_id) {
	return (bkt->value >> COUNT_LOC[type_id][count_index]) & ((1UL << COUNT_LEN[type_id][count_index]) - 1UL);
}

static inline __attribute__((always_inline)) void set_bucket_count(ec_bucket* bkt, int count_index, uint64_t count_value, const uint32_t type_id) {
	bkt->value = (bkt->value & ~(((1UL << COUNT_LEN[type_id][count_index]) - 1UL) << COUNT_LOC[type_id][count_index])) | ((count_value & ((1UL << COUNT_LEN[type_id][count_index]) - 1UL)) << COUNT_LOC[type_id][count_index]);
}


#define GET_HASH_VALUE_SENTENCE(key) int16_t final_key_len = 14;\
		h0_0 = (uint)(AwareHash((unsigned char *)key, 14, h[0], s[0], n[0]) % row_number);\
		h1_0 = (uint)(AwareHash((unsigned char *)key, 14, h[1], s[1], n[1]) % row_number); \
		h2_0 = (uint)(AwareHash((unsigned char *)key, 14, h[2], s[2], n[2]) % row_number);\
		h3_0 = (uint)(AwareHash((unsigned char *)key, 14, h[3], s[3], n[3]) % row_number); \
		h4_0 = (uint)(AwareHash((unsigned char *)key, 14, h[4], s[4], n[4]) % row_number); \
		uint8_t fp_0 = (SDBM((const unsigned char*)key, 14)  & FINGERPRINT_BITMASK);  \
		uint8_t fp_1 = (OAAT((const unsigned char*)key, 14)  & FINGERPRINT_BITMASK);  \
		uint8_t fp_2 = (BOB3((const unsigned char*)key, 14)  & FINGERPRINT_BITMASK);  \
		uint8_t fp_3 = (FNV32((const unsigned char*)key, 14)  & FINGERPRINT_BITMASK);  \
		uint8_t fp_4 = (BOB1((const unsigned char*)key, 14)  & FINGERPRINT_BITMASK);  \
		if (fp_0 == 0) fp_0 = 111 & FINGERPRINT_BITMASK;\
		if (fp_1 == 0) fp_1 = 127 & FINGERPRINT_BITMASK;\
		if (fp_2 == 0) fp_2 = 103 & FINGERPRINT_BITMASK;\
		if (fp_3 == 0) fp_3 = 85 & FINGERPRINT_BITMASK;\
		if (fp_4 == 0) fp_3 = 103 & FINGERPRINT_BITMASK;\
		h0_1 = ( h0_0 ^ (fp_0) ) % row_number; \
		h1_1 = ( h1_0 ^ (fp_1) ) % row_number; \
		h2_1 = ( h2_0 ^ (fp_2) ) % row_number; \
		h3_1 = ( h3_0 ^ (fp_3) ) % row_number; \
		h4_1 = ( h4_0 ^ (fp_4) ) % row_number; \
		uint hash[10] = {h0_0, h0_1, h1_0, h1_1, h2_0, h2_1, h3_0, h3_1, h4_0, h4_1}; \
		uint8_t fpSet[5] = {fp_0, fp_1, fp_2, fp_3, fp_4};



using namespace std;
template<int _bucket, int _row_number, int _cycle>
class Hermetis
{
public:
	uint cycle, bucket_num, row_number, maxloop,h1, h2, pos1, pos2;	//bucket_num indicates the number of buckets in each array
	uint h0_0, h0_1, h1_0, h1_1, h2_0, h2_1, h3_0, h3_1 ,h4_0, h4_1 ;
	int Delta_0 = 0, Delta_1 = 0, Delta_2 = 0, Delta_3 = 0 , Delta_4 = 0;
	double step;
	double Time_threshold;
	unsigned int clock_pos;
	unsigned int last_time;
	unsigned int NUM;
	Unit* counter;
	BOBHash *bobhash[5];		//Bob hash function
	// cm sketch
	#define CM_SKETCH_WIDTH 2
	#define CM_MAX_CNT 3
	uint cm_sketch_num;
	uint8_t *cm_sketch[CM_SKETCH_WIDTH];
	BOBHash *cm_hash[CM_SKETCH_WIDTH];
	std::map<uint64_t, uint64_t> distribution;

	uint64_t h[5], s[5], n[5];

	Hermetis() {
		bucket_num = _bucket;
		NUM = 0;
		row_number = calc_next_prime(_row_number);
		// row_number = _row_number;
		cycle = _cycle;
		step = (3.0 * row_number) / cycle;
		counter = new Unit [row_number];

		
		Time_threshold = (double)(cycle) / 3.0;
		clock_pos = 0;
        last_time = 0;

		for(int i = 0; i < row_number; i++){
			counter[i].count = new ec_bucket[bucket_num];
			memset(counter[i].count, 0, sizeof(ec_bucket) * bucket_num);
		}

		for (int i = 0; i < 5; i++) {
			bobhash[i] = new BOBHash(i + 1000);
		}

		int index = 0;
		for (int i = 0; i < 5; ++i) {
			h[i] = GenHashSeed(index++);
			s[i] = GenHashSeed(index++);
			n[i] = GenHashSeed(index++);
		}

		init_bucket_parameters();

	}
	
	// Insert the key into the counter
	inline int CM_insert(const char *key, const int16_t key_len = 0) {
		uint cm_idx[CM_SKETCH_WIDTH];
		for (int i = 0; i < CM_SKETCH_WIDTH; i++) {
			cm_idx[i] = (cm_hash[i]->run(key, key_len)) % cm_sketch_num;
			if ( cm_sketch[i][cm_idx[i]] != CM_MAX_CNT ) {
				cm_sketch[i][cm_idx[i]]++;
			}
		}
		return 0;
	}

	inline uint64_t CM_query(const char *key, const int16_t key_len = 0) {
		uint cm_idx[CM_SKETCH_WIDTH];
		int min_value = CM_MAX_CNT;
		for (int i = 0; i < CM_SKETCH_WIDTH; i++) {
			cm_idx[i] = (cm_hash[i]->run(key, key_len)) % cm_sketch_num;
			if (cm_sketch[i][cm_idx[i]] < min_value) {
				min_value = cm_sketch[i][cm_idx[i]];
			}
		}
		return (min_value < CM_MAX_CNT)? min_value: (100000);
	}

	// Copy a[0] --> b[0], a[1] --> b[1], ..., a[item_num-1] --> b[item_num-1]
	inline void copy_items_one_by_one(ec_bucket *dst, ec_bucket *src) {
		const uint32_t src_type_id = get_bucket_type_id(src);
		const uint32_t dst_type_id = get_bucket_type_id(dst);
		const uint32_t item_num = get_item_num_in_bucket_type(dst_type_id);
		for (int idx = 0; idx < item_num; idx++) {
			set_bucket_fingerprint(dst, idx, get_bucket_fingerprint(src, idx));
			set_bucket_time_tag(dst, idx, get_bucket_time_tag(src, idx, src_type_id), dst_type_id);
			set_bucket_count(dst, idx, get_bucket_count(src, idx, src_type_id), dst_type_id); 
			// Guarantee the copy is self
			if ( (1UL << COUNT_LEN[dst_type_id][idx]) - 1 < get_bucket_count(src, idx, src_type_id) ) {
				printf("The dst_type_id is %d, the idx is %d, the src_type_id is %d\n", dst_type_id, idx, src_type_id);
				assert(false);
			}
		}
	}

	// Copy a[1] --> b[0], a[2] --> b[1], ..., a[item_num] --> b[item_num-1]
	inline void copy_items_upflow(ec_bucket *dst, ec_bucket *src) {
		const uint32_t src_type_id = get_bucket_type_id(src);
		const uint32_t dst_type_id = get_bucket_type_id(dst);
		const uint32_t item_num = get_item_num_in_bucket_type(dst_type_id);
		for (int idx = 0; idx < item_num; idx++) {
			set_bucket_fingerprint(dst, idx, get_bucket_fingerprint(src, 1+idx));
			set_bucket_time_tag(dst, idx, get_bucket_time_tag(src, 1+idx, src_type_id), dst_type_id);
			set_bucket_count(dst, idx, get_bucket_count(src, 1+idx, src_type_id), dst_type_id);

			if ( COUNT_LEN[dst_type_id][idx] < COUNT_LEN[src_type_id][idx + 1] ) {
				printf("The dst_type_id is %d, the idx is %d, the src_type_id is %d\n", dst_type_id, idx, src_type_id);
				assert(false);
			}
		}
	}

	inline bool kick_to(int origin_hash_table_idx, uint32_t origin_slot_idx, uint8_t fingerprint_value, uint64_t count_value, uint8_t time_tag) {
	if (fingerprint_value == 0) {
			return true; // The kick out is NULL
		}
		ec_bucket *origin_bucket = counter[origin_slot_idx].count + (origin_hash_table_idx);
		uint32_t new_row_idx = ( origin_slot_idx ^ (fingerprint_value) ) % row_number;
		ec_bucket *dst_bucket = counter[new_row_idx].count + (origin_hash_table_idx);
	
		if (maxloop-- != 0) {
			uint32_t new_bucket_type_id = get_bucket_type_id(dst_bucket);
			uint8_t slot_num = get_item_num_in_bucket_type(new_bucket_type_id);
			uint min_slot_idx = 6;
			bool flag = false;
			for (int i = 0; i < slot_num; i++) {
				if (get_bucket_fingerprint(dst_bucket, i) == 0 &&  ( 1UL << COUNT_LEN[new_bucket_type_id][i]) > count_value) {
					set_bucket_count(dst_bucket, i, count_value, new_bucket_type_id);
					set_bucket_time_tag(dst_bucket, i, time_tag, new_bucket_type_id);
					set_bucket_fingerprint(dst_bucket, i, fingerprint_value);
					return true;
				}else if(( 1UL << COUNT_LEN[new_bucket_type_id][i]) > count_value && get_bucket_fingerprint(dst_bucket, i) != fingerprint_value){
					if(!flag ){
						min_slot_idx = i;
						flag = true;
					}
				}
			}
			if(!flag){
				return false; // The kicking fails
			}
			
			uint32_t new_type_id = get_bucket_type_id(dst_bucket);
			uint8_t new_time_tag = get_bucket_time_tag(dst_bucket, min_slot_idx, new_type_id);
			uint64_t new_count = get_bucket_count(dst_bucket, min_slot_idx, new_type_id);
			uint8_t new_fingerprint = get_bucket_fingerprint(dst_bucket, min_slot_idx);

			bool is_kick = kick_to(origin_hash_table_idx, new_row_idx, new_fingerprint, new_count,new_time_tag);
			if(is_kick){
				set_bucket_count(dst_bucket, min_slot_idx, count_value, new_bucket_type_id);
				set_bucket_time_tag(dst_bucket, min_slot_idx, time_tag, new_bucket_type_id);
				set_bucket_fingerprint(dst_bucket, min_slot_idx, fingerprint_value);
			}
			// std::cout <<"Hermetis kick_to result is is_kick ==  " << is_kick << std::endl;
			return is_kick;
		}
		
		 return false;// The kicking fails
	}
	
	inline bool solve_overflow_locally(ec_bucket* b, const int finger_idx, const uint32_t type_id, const uint32_t hash_table_idx, const uint32_t slot_idx, const uint32_t v = 1) { //overflow occur in entry_j in the bucket
		for ( int i = finger_idx + 1; i < get_item_num_in_bucket_type(type_id); i++ ) {
			// Try to find an empty entry with more bits
			if ( get_bucket_fingerprint(b, i) == 0 && (get_bucket_count(b, finger_idx, type_id) + v < (1UL << COUNT_LEN[type_id][i]))) {
				// Find an empty entry: move the data and do not change the type
				set_bucket_fingerprint(b, i, get_bucket_fingerprint(b, finger_idx));
				set_bucket_count(b, i, get_bucket_count(b, finger_idx, type_id) + v, type_id);
				set_bucket_time_tag(b, i, get_bucket_time_tag(b, finger_idx, type_id), type_id);
				
				set_bucket_fingerprint(b, finger_idx, 0);
				set_bucket_count(b, finger_idx, 0, type_id);
				set_bucket_time_tag(b, finger_idx, 0, type_id);
				return true;
			}
		}
		// Try to exchange with other slot in the same bucket
		uint64_t src_count = get_bucket_count(b, finger_idx, type_id) + v;
		uint8_t src_time_tag = get_bucket_time_tag(b, finger_idx, type_id);
		uint8_t src_fingerprint = get_bucket_fingerprint(b, finger_idx);
		for ( int dst_idx = finger_idx + 1; dst_idx < get_item_num_in_bucket_type(type_id); dst_idx++ ) {
			if ( get_bucket_count(b, dst_idx, type_id) < src_count - 1 ) {
				// Upper slot can be exchanged
				uint8_t dst_fingerprint = get_bucket_fingerprint(b, dst_idx);
				uint8_t dst_time_tag = get_bucket_time_tag(b, dst_idx, type_id);
				uint64_t dst_count = get_bucket_count(b, dst_idx, type_id);
				// Exchange
				set_bucket_fingerprint(b, dst_idx, src_fingerprint );
				set_bucket_time_tag(b, dst_idx, src_time_tag, type_id);
				set_bucket_count(b, dst_idx, src_count, type_id);
				set_bucket_fingerprint(b, finger_idx, dst_fingerprint);
				set_bucket_count(b, finger_idx, dst_count, type_id);
				set_bucket_time_tag(b, finger_idx, dst_time_tag, type_id);
				return true;
			}
		}
		// No empty entry found, we try to change it to other type of bucket
		ec_bucket new_bkt; new_bkt.value = 0;
		uint8_t least_finger, least_time_tag; uint64_t least_count;
		switch (type_id) {
			case 0:
			{
				bool flag_need_up = false;
				uint32_t next_type_id;
				if ( finger_idx != 0 ) {
					flag_need_up = true;
					next_type_id = (finger_idx == 3)? 1 : 2;
				}
				// It is the bottom type! We only kick out the least entry
				least_finger = get_bucket_fingerprint(b, 0);
				least_count = get_bucket_count(b, 0, type_id);
				least_time_tag = get_bucket_time_tag(b, 0, type_id);
				if (flag_need_up) {
					// Go right
					set_bucket_type_id(&new_bkt, next_type_id);
					copy_items_upflow(&new_bkt, b);
					b->value = new_bkt.value;
				} else {
					// Keep the same type if finger_idx == 0
					set_bucket_fingerprint(b, 0, 0);
					set_bucket_count(b, 0, 0, type_id); // Clear the least slot
					set_bucket_time_tag(b, 0, 0, type_id);
					// assert(least_count == 3 || least_count == 31);
				}
				if ( kick_to(hash_table_idx, slot_idx, least_finger,  least_count, least_time_tag) ) {
					return true;
				} else if ( !flag_need_up ) {
					// Keep the previous value
					set_bucket_fingerprint(b, 0, least_finger);
					set_bucket_count(b, 0, least_count, type_id);
					set_bucket_time_tag(b, 0, least_time_tag, type_id);
					return false;
				} else {
					#ifdef ENABLE_CM_SKETCH
					char cm_key[10];
					uint new_slot_idx = slot_idx;
					if (hash_table_idx == 1) {
						new_slot_idx = new_slot_idx ^ ((uint) least_finger);
					} 
					memset(cm_key, 0, sizeof(char) * 10);
					memcpy(cm_key, &new_slot_idx, sizeof(uint));
					memcpy(cm_key+sizeof(uint), &least_finger, sizeof(least_finger));
					uint64_t cm_value = CM_query(cm_key, sizeof(uint)+sizeof(char));
					#endif
				}
				break;
			}
			case 1 ... 5:
			{
				// Go down
				uint32_t possible_next_type_id = type_id + 1;
				set_bucket_type_id(&new_bkt, possible_next_type_id);
				uint32_t MAX_LEN_FOR_NEXT_TYPE = COUNT_LEN[possible_next_type_id][2];
				uint64_t MAX_COUNT_VALUE_FOR_NEXT = ((1ULL << MAX_LEN_FOR_NEXT_TYPE) - 2);
				if ( get_bucket_count(b, 2, type_id) <= MAX_COUNT_VALUE_FOR_NEXT ) {
					copy_items_one_by_one(&new_bkt, b);
					b->value = new_bkt.value;
					return true;
				} else {
					uint8_t out_finger = get_bucket_fingerprint(b, finger_idx);
					uint64_t out_count = get_bucket_count(b, finger_idx, type_id);
					uint8_t out_time_tag = get_bucket_time_tag(b, finger_idx, type_id);
					set_bucket_fingerprint(b, finger_idx, 0);
					set_bucket_count(b, finger_idx, 0, type_id);
					set_bucket_time_tag(b, finger_idx, 0, type_id);
					bool res = kick_to(hash_table_idx, slot_idx, out_finger, out_count,out_time_tag);
					if (res) {
						return true;
					} else {
						#ifdef ENABLE_CM_SKETCH
						char cm_key[10];
						uint new_slot_idx = slot_idx;
						if (hash_table_idx == 1) {
							new_slot_idx = new_slot_idx ^ ((uint) out_finger);
						} 
						memset(cm_key, 0, sizeof(char) * 10);
						memcpy(cm_key, &new_slot_idx, sizeof(uint));
						memcpy(cm_key+sizeof(uint), &out_finger, sizeof(out_finger));
						uint64_t cm_value = CM_query(cm_key, sizeof(uint)+sizeof(char));
						#endif
						return false;
					}
				}
				break;
			}
		
			default:
				assert(false);
		}
		return false;
	}
	
	inline bool aging_table_inversion(ec_bucket* b, const int finger_idx, const uint32_t type_id){
		if(type_id == 0) return true;
		// std::cout <<"Hermetis aging_table_inversion start ! type_id == " << type_id << std::endl;
		ec_bucket new_bkt; new_bkt.value = 0;
		set_bucket_fingerprint(b, finger_idx, 0);
		set_bucket_count(b, finger_idx, 0, type_id); // Clear the least slot
		set_bucket_time_tag(b, finger_idx, 0, type_id);
		switch (type_id){
			case 6:
			{
				uint32_t possible_inverse_type_id = type_id - 1;
				set_bucket_type_id(&new_bkt, possible_inverse_type_id);
				for (int i = 0,j = 0; i < 3; ){
					if(get_bucket_fingerprint(b, i) == 0){
						i++;
					}
					if(get_bucket_count(b, i, type_id) <= (1UL << COUNT_LEN[possible_inverse_type_id][j] - 1UL)){
						set_bucket_fingerprint(&new_bkt, j, get_bucket_fingerprint(b, i));
						set_bucket_count(&new_bkt, j, get_bucket_count(b,i,type_id), possible_inverse_type_id); 
						set_bucket_time_tag(&new_bkt, j, get_bucket_time_tag(b,i,type_id), possible_inverse_type_id); 
						i++,j++;
					}else{
						j++;
					}
				}
				b->value = new_bkt.value;
				break;
			}
			case 3 ... 5:{
				if(finger_idx){
					uint32_t possible_inverse_type_id = type_id - 1;
					set_bucket_type_id(&new_bkt, possible_inverse_type_id);
					for (int i = 2,j = 2; i >=0,j >=0; ){
						if(get_bucket_fingerprint(b, i) == 0){
							i--;
						}else{
							set_bucket_fingerprint(&new_bkt, j, get_bucket_fingerprint(b, i));
							set_bucket_count(&new_bkt, j, get_bucket_count(b,i,type_id), possible_inverse_type_id); 
							set_bucket_time_tag(&new_bkt, j, get_bucket_time_tag(b,i,type_id), possible_inverse_type_id); 
							i--,j--;
						}
					}
					b->value = new_bkt.value;
				}
				break;
			}
			case 2:{
				if(finger_idx == 1){
					uint32_t possible_inverse_type_id = 1;
					set_bucket_type_id(&new_bkt, possible_inverse_type_id);
					for (int i = 2,j = 2; i >=0,j >=0; ){
						if(get_bucket_fingerprint(b, i) == 0){
							i--;
						}else{
							set_bucket_fingerprint(&new_bkt, j, get_bucket_fingerprint(b, i));
							set_bucket_count(&new_bkt, j, get_bucket_count(b,i,type_id), possible_inverse_type_id); 
							set_bucket_time_tag(&new_bkt, j, get_bucket_time_tag(b,i,type_id), possible_inverse_type_id); 
							i--,j--;
						}
					}
					b->value = new_bkt.value;
				}else if(finger_idx == 2){
					uint32_t possible_inverse_type_id = 0;
					set_bucket_type_id(&new_bkt, possible_inverse_type_id);
					for (int i = 0,j = 0; i < 3; ){
						if(get_bucket_fingerprint(b, i) == 0){
							i++;
						}
						if(get_bucket_count(b, i, type_id) <= (1UL << COUNT_LEN[possible_inverse_type_id][j] - 1UL)){
							set_bucket_fingerprint(&new_bkt, j, get_bucket_fingerprint(b, i));
							set_bucket_count(&new_bkt, j, get_bucket_count(b,i,type_id), possible_inverse_type_id); 
							set_bucket_time_tag(&new_bkt, j, get_bucket_time_tag(b,i,type_id), possible_inverse_type_id); 
							i++,j++;
						}else{
							j++;
						}
					}
				}
				b->value = new_bkt.value;
				break;
			}
			case 1:{
				if(finger_idx == 2 || (get_bucket_count(b, 2, type_id)) <= (1UL << COUNT_LEN[0][3] - 1UL)){
					uint32_t possible_inverse_type_id = 0;
					set_bucket_type_id(&new_bkt, possible_inverse_type_id);
					for (int i = 0,j = 0; i < 3; ){
						if(get_bucket_fingerprint(b, i) == 0){
							i++;
						}
						if(get_bucket_count(b, i, type_id) <= (1UL << COUNT_LEN[possible_inverse_type_id][j] - 1UL)){
							set_bucket_fingerprint(&new_bkt, j, get_bucket_fingerprint(b, i));
							set_bucket_count(&new_bkt, j, get_bucket_count(b,i,type_id), possible_inverse_type_id); 
							set_bucket_time_tag(&new_bkt, j, get_bucket_time_tag(b,i,type_id), possible_inverse_type_id); 
							i++,j++;
						}else{
							j++;
						}
					}
					b->value = new_bkt.value;
				}
				break;
			}
		}
		return true;
	}

	inline bool plus(ec_bucket* b, const int finger_idx, const uint32_t type_id, const int hash_table_idx, const uint32_t slot_idx,  const int cur_timetag,  const uint32_t v = 1) { //try to plus entry_j in the bucket, return true if no overflow happens
		uint64_t original_val = get_bucket_count(b, finger_idx, type_id);
		int original_time_tag = get_bucket_time_tag(b, finger_idx, type_id);
		if((cur_timetag - original_time_tag + 8) % 4 == 3){
			set_bucket_time_tag(b, finger_idx, cur_timetag, type_id);
			set_bucket_count(b, finger_idx, v, type_id);
		}else{
			const uint64_t max_cnt_val = (1UL << COUNT_LEN[type_id][finger_idx]);
			if ( v + original_val >= max_cnt_val ) { 
				bool res = solve_overflow_locally(b, finger_idx, type_id, hash_table_idx, slot_idx, v);
				return res;
			}
			original_val += v;
			set_bucket_count(b, finger_idx, original_val, type_id);
		
		}
		return true;
	}
	
	inline void flush(unsigned int num){
		int cur_time_tag_0 = ((int)((num + Delta_0) / Time_threshold)) % 4;
		int cur_time_tag_1 = ((int)((num + Delta_1) / Time_threshold)) % 4;
		int cur_time_tag_2 = ((int)((num + Delta_2) / Time_threshold)) % 4;
		int cur_time_tag_3 = ((int)((num + Delta_3) / Time_threshold)) % 4;
		int cur_time_tag_4 = ((int)((num + Delta_4) / Time_threshold)) % 4;
		int timg_tag[5] = {cur_time_tag_0, cur_time_tag_1, cur_time_tag_2, cur_time_tag_3,cur_time_tag_4};

		for(int slice = 0;slice < 4;++slice){
			for(int Row = 0;Row < row_number;++Row){
				int cur_timetag = timg_tag[slice];
				ec_bucket *b = counter[Row].count + (slice);
				uint32_t type_id = get_bucket_type_id(b);
				uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
				for (int j = fingerprint_num-1; j >= 0; j--) {
				// for (int j = 0; j < fingerprint_num; j++) {
					if (get_bucket_fingerprint(b, j)  && (cur_timetag - (int)get_bucket_time_tag(b, j, type_id)) % 4 == 3) {
						// std::cout << "Aging clear happens!  " << get_bucket_count(b, j, type_id) << std::endl;
						set_bucket_fingerprint(b, j, 0);
						set_bucket_count(b, j, 0, type_id);
						set_bucket_time_tag(b, j, 0, type_id);
						if(type_id != 0) aging_table_inversion(b,j,type_id);
					} 
				}
			}
		}

	}

	inline void Clock_Go(unsigned int num){
		int cur_time_tag_0 = ((int)((num + Delta_0) / Time_threshold)) % 4;
		int cur_time_tag_1 = ((int)((num + Delta_1) / Time_threshold)) % 4;
		int cur_time_tag_2 = ((int)((num + Delta_2) / Time_threshold)) % 4;
		int cur_time_tag_3 = ((int)((num + Delta_3) / Time_threshold)) % 4;
		int cur_time_tag_4 = ((int)((num + Delta_4) / Time_threshold)) % 4;
		int timg_tag[5] = {cur_time_tag_0, cur_time_tag_1, cur_time_tag_2, cur_time_tag_3,cur_time_tag_4};

		for(int slice = 0;slice < 5;++slice){
			int cur_timetag = timg_tag[slice];
			for(auto itr = hotKey_Pos_st.begin();itr != hotKey_Pos_st.end();++itr){
				ec_bucket *b = counter[*itr].count + (slice);
				uint32_t type_id = get_bucket_type_id(b);
				uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
				for (int j = fingerprint_num-1; j >= 0; j--) {
				// for (int j = 0; j < fingerprint_num; j++) {
					if (get_bucket_fingerprint(b, j)  && (cur_timetag - (int)get_bucket_time_tag(b, j, type_id)) % 4 == 3) {
						// std::cout << "*itr == " << *itr << std::endl;
						set_bucket_fingerprint(b, j, 0);
						set_bucket_count(b, j, 0, type_id);
						set_bucket_time_tag(b, j, 0, type_id);
						if(type_id != 0) aging_table_inversion(b,j,type_id);
					} 
				}
			}
		}	

		for(;last_time < num * step + 1;++last_time){
			for(int slice = 0;slice < 4;++slice){
				int cur_timetag = timg_tag[slice];
				ec_bucket *b = counter[clock_pos].count + (slice);
				uint32_t type_id = get_bucket_type_id(b);
				uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
				for (int j = fingerprint_num-1; j >= 0; j--) {
				// for (int j = 0; j < fingerprint_num; j++) {
					if (get_bucket_fingerprint(b, j)  && (cur_timetag - (int)get_bucket_time_tag(b, j, type_id)) % 4 == 3) {
						// std::cout << "Aging clear happens!  " << get_bucket_count(b, j, type_id) << std::endl;
						set_bucket_fingerprint(b, j, 0);
						set_bucket_count(b, j, 0, type_id);
						set_bucket_time_tag(b, j, 0, type_id);
						if(type_id != 0) aging_table_inversion(b,j,type_id);
					} 
				}
			}
			clock_pos = (clock_pos + 1) % row_number;
		}
		
		
	}

	void Insert(const char *key,  unsigned int num, const uint32_t v = 1, const int16_t key_len = 14) {
		NUM = num;
		int cur_time_tag_0 = ((int)((num + Delta_0) / Time_threshold)) % 4;
		int cur_time_tag_1 = ((int)((num + Delta_1) / Time_threshold)) % 4;
		int cur_time_tag_2 = ((int)((num + Delta_2) / Time_threshold)) % 4;
		int cur_time_tag_3 = ((int)((num + Delta_3) / Time_threshold)) % 4;
		int cur_time_tag_4 = ((int)((num + Delta_4) / Time_threshold)) % 4;
		int timg_tag[5] = {cur_time_tag_0, cur_time_tag_1, cur_time_tag_2, cur_time_tag_3, cur_time_tag_4};
		GET_HASH_VALUE_SENTENCE(key);
		bool flag = false;
		int empty_jj, empty_type_id; 
		ec_bucket* empty_bucket;
		for(int hash_table_idx = 0;hash_table_idx < 5;++hash_table_idx){
			maxloop = MAXLOOP;
			flag = false;
			int cur_timetag = timg_tag[hash_table_idx];
			uint8_t fp = fpSet[hash_table_idx];
			for (int i = 0; i < 2; i++) {
				ec_bucket *b = counter[hash[2 * hash_table_idx + i]].count + (hash_table_idx);
				uint32_t type_id = get_bucket_type_id(b);
				uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
				for (int j = fingerprint_num-1; j >= 0; j--) {
					if ( get_bucket_fingerprint(b, j) == fp ) {
						if((cur_timetag  + 4 - (int)get_bucket_time_tag(b, j, type_id)) % 4 == 3){
							set_bucket_count(b, j, v, type_id);
							set_bucket_time_tag(b, j, cur_timetag, type_id);
							break;
						}else{
							bool res = plus(b, j, type_id, hash_table_idx,  hash[2 * hash_table_idx + i], cur_timetag, v);
							break;
						}
						if(get_bucket_count(b, j, type_id) > 100){
							hotKey_Pos_st.emplace(hash[2 * hash_table_idx + i]);
						}
						
					} else if ( !flag && get_bucket_fingerprint(b, j) == 0) {
						empty_bucket = b;
						empty_jj = j;
						empty_type_id = type_id;
						flag = true; // We have an empty bucket
					}
					}
			}
			if (flag) {
				// Insert the key into the empty bucket
				set_bucket_fingerprint(empty_bucket, empty_jj, fp);
				set_bucket_count(empty_bucket, empty_jj, v, empty_type_id);
				set_bucket_time_tag(empty_bucket, empty_jj, cur_timetag, empty_type_id);
			} 
		}
		return ;
	}

	double zero(){
		double cnt[2]={0.0, 0.0};
		for (int i=0; i<2; i++) {
			for (int j=0; j<bucket_num; j++) {
				ec_bucket *b = counter[i].count + j;
				const uint64_t type_id = get_bucket_type_id(b);
				const uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
				int flag = 1;
				for (int k=0; k< fingerprint_num; k++) {
					if ( get_bucket_fingerprint(b, k) != 0 ) { flag = 0; break; }
				}
				if (1 == flag) cnt[i]+=1; // The bucket is empty
			}
		}
		return (bucket_num-cnt[0])/bucket_num * (bucket_num-cnt[1])/bucket_num;
	}	
	
	int Query(const char *key, unsigned int num, const int16_t key_len = 0,bool debug = false) {
		GET_HASH_VALUE_SENTENCE(key);
		int cur_time_tag_0 = ((int)((num + Delta_0) / Time_threshold)) % 4;
		int cur_time_tag_1 = ((int)((num + Delta_1) / Time_threshold)) % 4;
		int cur_time_tag_2 = ((int)((num + Delta_2) / Time_threshold)) % 4;
		int cur_time_tag_3 = ((int)((num + Delta_3) / Time_threshold)) % 4;
		int cur_time_tag_4 = ((int)((num + Delta_4) / Time_threshold)) % 4;
		int timg_tag[5] = {cur_time_tag_0, cur_time_tag_1, cur_time_tag_2, cur_time_tag_3, cur_time_tag_4};
		bool flag = false;
		int result = 20000;
		vector<int> resultSet;
		uint64_t min_value = UINT64_MAX; uint64_t table_min[2], max_value = 0;
		__builtin_prefetch(counter[hash[0]].count + 0, 0, 2);
		__builtin_prefetch(counter[hash[1]].count + 1, 0, 2);
		__builtin_prefetch(counter[hash[2]].count + 0, 0, 2);
		__builtin_prefetch(counter[hash[3]].count + 1, 0, 2);
		__builtin_prefetch(counter[hash[4]].count + 0, 0, 2);
		__builtin_prefetch(counter[hash[5]].count + 1, 0, 2);
		__builtin_prefetch(counter[hash[6]].count + 0, 0, 2);
		__builtin_prefetch(counter[hash[7]].count + 1, 0, 2);
		__builtin_prefetch(counter[hash[8]].count + 0, 0, 2);
		__builtin_prefetch(counter[hash[9]].count + 1, 0, 2);
		
		resultSet.clear();

		for(int slice = 0;slice < 5;++slice){
			flag = false;
			int cur_timetag = timg_tag[slice];
			uint8_t fp = fpSet[slice];
			for (int i = 0; i < 2; i++) {
				ec_bucket *b = counter[hash[2 * slice + i]].count + (slice);
				uint32_t type_id = get_bucket_type_id(b);
				uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
				for (uint8_t fpt_idx = 0; fpt_idx < fingerprint_num; fpt_idx++) {
					const uint8_t stored_fingerprint = get_bucket_fingerprint(b, fpt_idx);
					const uint64_t stored_count = get_bucket_count(b, fpt_idx, type_id);
					if ( stored_fingerprint == fp ) {
						result = min(result,(int)stored_count);
						resultSet.push_back((int)stored_count);

					}
					if (!flag && stored_fingerprint == 0) {
						flag = true;
					}
					if (flag) {continue;}  
					if (stored_fingerprint != 0 && min_value > stored_count) { min_value = stored_count;}
				}
			}

			std::sort(resultSet.begin(), resultSet.end());
			int max_pos = resultSet.size() - 1;
			if(resultSet.size() == 0) return 0;
			if(resultSet[max_pos] - resultSet[0] < 4){
				return resultSet[max_pos] ;
			}else{
				if(resultSet.size() & 1){
					return resultSet[resultSet.size() >> 1];
				}else{
					int res = (resultSet[resultSet.size() >> 1] + resultSet[resultSet.size() >> 1] - 1) >> 1;
					if(res == 0) return resultSet[resultSet.size() - 1] ;
					return res;
				}
			}
			
			if(result != 20000) return result;
			if (flag) { return 0; } 
			return min_value;
		}
	return 0;
	}
	
	//memeory access
	int Mem(const char *key, const int16_t key_len = 0) {
		GET_HASH_VALUE_SENTENCE(key);
		for (uint8_t i = 0; i < 2; i++) {
			uint8_t fp = fpSet[0];
			ec_bucket* b = counter[hash[i]].count + hash[2+i];
			const uint64_t type_id = get_bucket_type_id(b);
			const uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
			for (uint8_t j = 0; j < fingerprint_num; j++) {
				const uint8_t stored_fingerprint = get_bucket_fingerprint(b, j);
				if (stored_fingerprint == fp) { return 1;}
			}
		}
		return 2;
	}

	// the use ratio of the bucket
	double Ratio() {
		int used_num = 0;
		int total_slot = 0;
		unordered_map<uint16_t, uint64_t> type_count;
		type_count.clear();
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < bucket_num; j++) {
				ec_bucket* b = counter[i].count + j;
				const uint32_t type_id = get_bucket_type_id(b);

				if (type_count.find(type_id) == type_count.end()) {
					type_count[type_id] = 1;
				} else {
					type_count[type_id]++;
				}

				const uint32_t slot_num = get_item_num_in_bucket_type(type_id);
				total_slot += slot_num;
				const uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
				for (int k = 0; k < fingerprint_num; k++) {
					if ( get_bucket_count(b, k, type_id) != 0 ) { used_num++; }
				}
			}
		}

		printf("The bucket number is %d\n", 2*bucket_num);

		for (int i = 0;  i < 12; i++) {
			if (type_count.find(i) != type_count.end()) {
				printf("type %d: %d with ratio %.2f\%\n", i, type_count[i], ((double) type_count[i] ) * 100 / 2.0 / bucket_num);
			}
		}

		return used_num / (total_slot * 1.0);
	}

	void dump_to_file(FILE* fp) {
		// Write the result to fp
		unordered_map<uint16_t, uint64_t> type_count;
		type_count.clear();
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < bucket_num; j++) {
				ec_bucket* b = counter[i].count + j;
				const uint32_t type_id = get_bucket_type_id(b);
				const uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
				// Group_id, bucket_addr, bucket_type, fingerprint_num, fingerprint_1, count_1, fingerprint_2, count_2, ..., fingerprint_n, count_n
				fprintf(fp, "%d\t%d\t%d\t%d\t", i, j, type_id, fingerprint_num);
				for (int k = 0; k < fingerprint_num; k++) {
					uint32_t fingerprint = get_bucket_fingerprint(b, k);
					uint64_t cnt = get_bucket_count(b, k, type_id);
					// write to the file
					fprintf(fp, "%d\t%d\t", fingerprint, cnt);
				}
				fprintf(fp, "\n");
			}
		}
	}


	//delete the bucket
	void Delete(char *key, const int16_t key_len = 0) {
		GET_HASH_VALUE_SENTENCE(key);
		for (int i = 0; i < 2; i++) {
			uint8_t fp = fpSet[0];
			ec_bucket* b = counter[hash[i]].count + hash[2+i];
			const uint32_t type_id = get_bucket_type_id(b);
			const uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
			for (int j = 0; j < fingerprint_num; j++) {
				const uint8_t stored_fingerprint = get_bucket_fingerprint(b, j);
				if (stored_fingerprint == fp) {
					uint64_t origin_value = get_bucket_count(b, j, type_id);
					set_bucket_count(b, j, origin_value - 1, type_id);
					if ( origin_value == 1 ) {
						set_bucket_fingerprint(b, j, 0);
					}
					return;
				}
			}
		}
	}


	uint32_t get_distribution() {
		uint32_t nonEmpeyNum = 0;
		int cur_time_tag_0 = ((int)((NUM + Delta_0) / Time_threshold)) % 4;

		for (int j=0; j<row_number; j++) {
			bool flag = false;
				ec_bucket* b = counter[j].count + 0;
				const uint64_t type_id = get_bucket_type_id(b);
				const uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
				for (int k=0; k < fingerprint_num; k++) {
					if ( !flag && get_bucket_count(b, k, type_id) > 0 && get_bucket_fingerprint(b, k) != 0 && (cur_time_tag_0  + 4 - (int)get_bucket_time_tag(b, k, type_id)) % 4 != 3 ) {
						flag = true; 
					}
				}
				if(flag){
					nonEmpeyNum++;
					}
		}

        distribution.clear();
		for (int j=0; j<row_number; j++) {
			ec_bucket* b = counter[j].count + 0;
			const uint64_t type_id = get_bucket_type_id(b);
			const uint8_t fingerprint_num = get_item_num_in_bucket_type(type_id);
			for (int k=0; k < fingerprint_num; k++) {
				if ( get_bucket_count(b, k, type_id) != 0) {
					const uint64_t value = get_bucket_count(b, k, type_id);
					if (value != 0 && distribution.find(value) == distribution.end())
						distribution.emplace(value, 1);
					else  
						distribution[value] += 1;
						// cout <<"distribution[value] is " << distribution[value] << endl;
					}
				}
		}


		return nonEmpeyNum;
	}

	uint32_t get_cardinality() {
		uint32_t nonEmpeyNum = get_distribution();
		int cardinality = 0;
		double counter_num;
		for (auto iter = distribution.begin(); iter != distribution.end(); ++iter) {
			if(iter->first != 0){
				cardinality += iter->second;
			}
		}

		
		double rho = (double)cardinality / nonEmpeyNum;
		double rate = ((double)row_number - nonEmpeyNum) / row_number;

		distribution.clear();
		return row_number * log(1 / rate) ;
	}

	std::pair<uint64_t, double> get_entropy() {
		uint64_t total = 0;
		double entropy = 0;
		get_distribution();
		for (auto iter = distribution.begin(); iter != distribution.end(); ++iter) {
			if(iter->first != 0){
				total += iter->first * iter->second;
				entropy += iter->first * iter->second * log2(iter->first);
			}
		}
		return std::make_pair(total, entropy);
	}

	size_t get_memory_usage() {
		return 8 * 5 * row_number;
	}

	~Hermetis() {

	}
	};
#endif//_Hermetis_H
