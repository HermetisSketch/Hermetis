#ifndef _COMMON_H
#define _COMMON_H

#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>
#include <array>
#include <chrono>

#include "AwareHash.h"
#include "parameters.h"

typedef const char* Key_t;
typedef std::array<uint32_t, BUCKET_NUM> Hist_t;

typedef struct five_tuple {
	char str[CHARKEY_LEN];
	bool operator < (const five_tuple& a) const {
		bool flag = false;
		for (int i = 0; i < CHARKEY_LEN; ++i) {
			if (this->str[i] == a.str[i])
				continue;
			else if (this->str[i] < a.str[i]) {
				flag = 1;
				break;
			}
			else
				return false;
		}
		return flag;
	}
	bool operator == (const five_tuple& a) const {
		for (int i = 0; i < CHARKEY_LEN; ++i) {
			if (this->str[i] == a.str[i])
				continue;
			else
				return false;
		}
		return true;
	}
	bool isEmpty() {
		for (int i = 0; i < CHARKEY_LEN; ++i) {
			if (str[i] == 0)
				continue;
			else
				return false;
		}
		return true;
	}
	five_tuple() {
		memset(str, 0, sizeof(str));
	}
	five_tuple(char * s) {
		memcpy(str, s, CHARKEY_LEN);
	}
	five_tuple(Key_t s) {
		memcpy(str, s, CHARKEY_LEN);
	}
} five_tuple;

typedef struct data_t{
	five_tuple key;		//the five tuple
	uint16_t length;	//IP packet length
	unsigned int timestamp;
	
	data_t(char *s) {
		key = five_tuple(s);
		length = *(uint16_t *)(s + CHARKEY_LEN);
		timestamp = 0;
	}
	data_t(Key_t s) {
		key = five_tuple(s);
		length = *(uint16_t *)(s + CHARKEY_LEN);
		timestamp = 0;
	}
	data_t() {
		length = 0;
		timestamp = 0;
	}
} data_t;


typedef struct new_data_t{
	uint64_t key;		//the five tuple
	uint8_t bid;		//bucket id
	
	new_data_t(uint64_t k, uint8_t bucket_id) {
		key = k;
		bid = bucket_id;
	}
} new_data_t;

#endif


