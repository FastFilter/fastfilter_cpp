/*
Copyright (c) 2019 Advanced Micro Devices, Inc.
 
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Author: Alex D. Breslow 
        Advanced Micro Devices, Inc.
        AMD Research
*/
#ifndef _VECTOR_TYPES_H
#define _VECTOR_TYPES_H

#include <array>

typedef __uint128_t uint128_t;
//typedef uint64_t atom_t; // Tested up to __uint128_t inclusive
//typedef uint64_t hash_t;
typedef uint8_t counter_t; // Used only in one implementation of scans
//typedef uint64_t keys_t; // IMPORTANT: C++ has a key_t in its implementation
                         // of <random>.  Make sure you use this one.
#define CCF_KEY_SIZE 8

#if CCF_KEY_SIZE == 4
  typedef uint32_t atom_t;
  typedef uint32_t hash_t;
  typedef uint32_t keys_t;
  constexpr uint64_t _N = 8;
#elif CCF_KEY_SIZE == 8
  typedef uint64_t atom_t;
  typedef uint64_t hash_t;
  typedef uint64_t keys_t;
  constexpr uint64_t _N = 4;
#else
  #error "Only CCF_KEY_SIZE 4 and 8 are currently supported"
#endif

// Could add the bucket id and fingerprint since they are associated with this
struct StoreParams{ // TODO: Optimize this a bit
  hash_t block_id;              // Could go
  uint16_t counter_index;       // Could go
  counter_t bucket_start_index; // Stay
  counter_t elements_in_block;  // Stay
  counter_t counter_value;         // Could go
};


constexpr uint_fast64_t batch_size = 128;
typedef std::array<atom_t, batch_size> ar_atom;
typedef std::array<uint8_t, batch_size> ar_u8;
typedef std::array<uint16_t, batch_size> ar_u16;
typedef std::array<uint32_t, batch_size> ar_u32;
typedef std::array<counter_t, batch_size> ar_counter;
typedef std::array<hash_t, batch_size> ar_hash;
typedef std::array<keys_t, batch_size> ar_key;
typedef std::array<StoreParams, batch_size> ar_store_params;

struct StoreParamsSOA{
  ar_hash block_ids;
  ar_u16 counter_indexes;
  ar_counter bucket_start_indexes;
  ar_counter elements_in_blocks;
  ar_counter counter_values;
};


// TODO: Check that this is GNU source
typedef uint64_t vN_u64 __attribute__ ((vector_size(sizeof(uint64_t) * _N)));
typedef uint8_t vN_u8 __attribute__ ((vector_size (sizeof(uint8_t) * _N)));
typedef uint16_t vN_u16 __attribute__ ((vector_size (sizeof(uint16_t) * _N)));
typedef uint32_t vN_u32 __attribute__ ((vector_size (sizeof(uint32_t) * _N)));
typedef atom_t vN_atom __attribute__ ((vector_size (sizeof(atom_t) * _N)));
typedef keys_t vN_key __attribute__((vector_size (sizeof(keys_t) * _N)));
typedef hash_t vN_hash __attribute__((vector_size (sizeof(hash_t) * _N)));

// For some reason using atom_t for counters in the vectorized case is faster.
typedef atom_t vN_counter __attribute__ ((vector_size (sizeof(atom_t) * _N)));

#endif
