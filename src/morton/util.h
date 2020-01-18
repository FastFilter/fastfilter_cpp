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
#ifndef _UTIL_H
#define _UTIL_H

#include <cstdint>
#include <string>
#include <sstream>
#include <cmath>

#include <iostream>

#include "vector_types.h"

// FIXME: Put guards around this
// For BMI2 pdep instruction
#ifdef __BMI2__
#include "x86intrin.h"
#endif

namespace util{

  template<class INT_TYPE>
  inline std::string bin_string(INT_TYPE integer, uint32_t spacing){
    std::stringstream ss;
    for(int32_t i = sizeof(integer) * 8 - 1; i > -1; i--){
      ss << ((integer >> i) & 1 ? '1' : '0');
      if(i % spacing == 0){
        ss << ' '; 
      }
    }
    return ss.str();
  }

  template<class INT_TYPE>
  inline std::string bin_string(INT_TYPE integer){
    return bin_string<INT_TYPE>(integer, 8 * sizeof(INT_TYPE));
  }

  // This could be implemented using fancy binary arithmatic or builtins, 
  // but this probably suffices if the integer is known at compile time.
  constexpr inline uint32_t log2ceil(uint32_t integer){
    //return ceil(log2(integer));
    return 32u - __builtin_clz(integer - 1u);
  }

  // See https://lemire.me/blog/2016/06/27
  // These functions implement a fast alternative to the modulo reduction.
  // The algorithm is presented by Professor Daniel Lemire of the University
  // of Quebec in his outstanding blog, which is under a Creative Commons 
  // Attribution-ShareAlike 3.0 Unported License. 
  // See https://creativecommons.org/licenses/by-sa/3.0/us/ and
  // https://lemire.me/blog/terms-of-use/.
  template<typename T>
  inline T fast_mod_alternative(T raw_hash, T modulus, T hash_width_in_bits);

  template<>
  inline uint64_t fast_mod_alternative<uint64_t>(uint64_t raw_hash, 
    uint64_t modulus, uint64_t hash_width_in_bits){
    return (static_cast<__uint128_t>(raw_hash) * modulus) >> hash_width_in_bits;
  }

  template<>
  inline uint32_t fast_mod_alternative<uint32_t>(uint32_t raw_hash, 
    uint32_t modulus, uint32_t hash_width_in_bits){
    return (static_cast<__uint64_t>(raw_hash) * modulus) >> hash_width_in_bits;
  }

  template<class TN, class T> 
  inline TN fast_mod_alternativeN(TN raw_hashes, T modulus);

	template<class ARRAY_TYPE>
	inline void print_array(const std::string& name, const ARRAY_TYPE& array){
		std::cout << name << " [ ";
    for(uint32_t i = 0; i < batch_size; i++){
      std::cout << static_cast<uint32_t>(array[i]) << " ";
		}
		std::cout << "]\n";
	}

  template<>
  inline vN_u32 fast_mod_alternativeN<vN_u32, uint32_t>(vN_u32 raw_hashes, uint32_t modulus){
    for(uint32_t i = 0; i < _N; i++){
      static_assert(_N <= 8, "Vector width exceeds AVX/AVX2's 256-bit vector width\n");
      raw_hashes[i] = static_cast<uint32_t>((static_cast<__uint64_t>(raw_hashes[i]) * modulus) >> 32U);
    }
    return raw_hashes;
  }

} // End of util namespace


// FIXME: Not yet tested
std::ostream& operator<<(std::ostream& os, __uint128_t integer){
  std::stringstream ss;
  __uint128_t sqrt_power10 = static_cast<__uint128_t>(10000000000000000000ull);
  // log10 of 2^127 is between 38 and 39, so start with 38 zeros
  __int128_t power10 = sqrt_power10 * sqrt_power10;
  while(static_cast<__uint128_t>(0) / power10 == 0){
    power10 /= 10;
  }
  while(power10 != 0){
    uint32_t digit = integer / power10;
    os << static_cast<uint32_t>(digit);
    integer -= power10 * (digit);
    power10 /= 10;
  }
  return os;
}


#endif
