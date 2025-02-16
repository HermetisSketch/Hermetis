#ifndef DATA_H
#define DATA_H

#include "definition.h"
#include "hash_class.h"
#include "parameters.h"

class Data{
public:
    unsigned char str[DATA_LEN];
    int timestamp;
    Data& operator = (Data an);
};

bool operator < (Data bn, Data an);
bool operator == (Data bn, Data an);

// class My_Hash{
// public:
//     size_t operator()(const Webget_Tuple dat) const{
//         uint64_t bid = dat.length / BUCKET_WIDTH;
//         new_data_t tmp_key((Key_t)(dat.key).str, bid);
//         // return RSHash(((const char*))tmp_key.str, 14);
//         return RSHash((const unsigned char*)tmp_key.str, 14);
//     }
// };

#endif // DATA_H
