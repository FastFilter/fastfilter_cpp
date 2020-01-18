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

// Some functions for hashing

#ifndef _HASH_UTIL_H
#define _HASH_UTIL_H

#include <random>

#include "vector_types.h"
#include "test_util.h"

struct BitMixMurmur {
  inline hash_t operator()(keys_t key) const{
    #if CCF_KEY_SIZE == 4
      return hash32N(key);
    #else
      return hash64N(key);
    #endif
  } 

  // TODO: Implement with template specialization
  template<class T> 
  inline T hashN(T ks) const{
   #if CCF_KEY_SIZE == 4
     return hash32N(ks);
   #else
     return hash64N(ks);
   #endif
  }

  // Based on code in the public domain (MurmurHash3a)
  // See https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp 
  // 92cf370
  // The disavowing of the copyright is reproduced below and applies only to MurmurHash3:
  //-----------------------------------------------------------------------------
  // MurmurHash3 was written by Austin Appleby, and is placed in the public
  // domain. The author hereby disclaims copyright to this source code.
  // Template it to make it usable with vector types such as vN_u64
  template <class T>
  inline T hash64N(T ks) const{ // Bit mix from MurmurHash64/CLHash
    ks ^= ks >> 33;
    ks *= 0xff51afd7ed558ccdULL;
    ks ^= ks >> 33;
    ks *= 0xc4ceb9fe1a85ec53ULL;
    ks ^= ks >> 33;
    return ks;
  }
  // Based on code in the public domain (MurmurHash3a)
  template <class T>
  inline T hash32N(T ks) const{ // Bit mix from MurmurHash32
    ks ^= ks >> 16;
    ks *= 0x85ebca6b;
    ks ^= ks >> 13;
    ks *= 0xc2b2ae35;
    ks ^= ks >> 16;
    return ks;
  }
};

#endif
