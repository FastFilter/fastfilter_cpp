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
#ifndef _BF_H
#define _BF_H

// Author: Alex Breslow
// Description: A simple specialized implementation of a cache-blocked
// Bloom filter.  See "Cache-, Hash- and Space-Efficient Bloom Filters" by
// Putze et al. in JEA'09
// URL: https://dl.acm.org/citation.cfm?id=1594230

#include <array>

#define LOG2(x) (x < 16 ? 3 : x < 32 ? 4 : x < 64 ? 5 : x < 128 ? 6 : x < 256 ? 7 : x < 512 ? 8 : -1)

namespace BlockedBF{
  using slot_type = uint32_t;
  struct Bucket{
    slot_type f1, f2, f3, f4;
  };
  template<uint64_t T_NUM_BUCKETS>
  struct BloomFilter{
    static_assert(__builtin_popcountll(T_NUM_BUCKETS) == 1, "BloomFilter must "
      "be a power of 2 buckets in size");
    std::array<Bucket, T_NUM_BUCKETS> buckets{};

    // Check if item is in the filter, then insert it.  Return if you found it.
    inline bool contains_and_update(const hash_t item){
      constexpr slot_type one = 1;
      constexpr slot_type slot_width_bits = sizeof(slot_type) * 8;
      constexpr slot_type log2_slot_width = LOG2(slot_width_bits);
      constexpr slot_type shift0 = 0;
      constexpr slot_type shift1 = log2_slot_width;
      constexpr slot_type shift2 = log2_slot_width * 2;
      constexpr slot_type shift3 = log2_slot_width * 3;
      constexpr slot_type bucket_shift = log2_slot_width * 4;
      // Multiply by MurmurHash constant
      uint64_t hashes = 0xff51afd7ed558ccdULL * item;
      slot_type h1 = (hashes >> shift0) % slot_width_bits;
      slot_type h2 = (hashes >> shift1) % slot_width_bits;
      slot_type h3 = (hashes >> shift2) % slot_width_bits;
      slot_type h4 = (hashes >> shift3) % slot_width_bits;
      slot_type bucket_id = (hashes >> bucket_shift) % (T_NUM_BUCKETS);
      bool conflict_present = (buckets[bucket_id].f1 >> h1) &
        (buckets[bucket_id].f2 >> h2) &
        (buckets[bucket_id].f3 >> h3) &
        (buckets[bucket_id].f4 >> h4) & one;
      buckets[bucket_id].f1 |= one << h1;
      buckets[bucket_id].f2 |= one << h2;
      buckets[bucket_id].f3 |= one << h3;
      buckets[bucket_id].f4 |= one << h4;
      return conflict_present;
    }
  };
}

#endif
