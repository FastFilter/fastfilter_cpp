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
#ifndef _COMPRESSED_CUCKOO_CONFIG_H
#define _COMPRESSED_CUCKOO_CONFIG_H

namespace CompressedCuckoo{
  // See vector_types.h for more types and for tuning the vector width and
  // atom types
  const bool g_cache_aligned_allocate = true;
  const size_t g_cache_line_size_bytes = 64;  // Change this as necessary
  const uint64_t stash_prefix_tag_len = 4;
  
  // Allows for up to 255 items per block
  const uint8_t max_fullness_counter_width = 8;
  constexpr atom_t one = static_cast<atom_t>(1);
  
  enum struct AlternateBucketSelectionMethodEnum{
    TABLE_BASED_OFFSET,
    FUNCTION_BASED_OFFSET,
    FAN_ET_AL_PARTIAL_KEY // Only use this if you can guarantee the total buckets 
                          // in the filter is a power of two
  };

  enum struct InsertionMethodEnum{
    FIRST_FIT,
    TWO_CHOICE,
    HYBRID_SIMPLE,
    HYBRID_PIECEWISE, // Starts off as first-fit and then transitions to two choice
                  // once you hit a certain load factor 
    FIRST_FIT_OPT, // Transitions between two implementations of first-fit
  };
 
  enum struct CounterReadMethodEnum{
    READ_SIMPLE,
    READ_CROSS,
    READ_RAW,  // If counters are always in atom 0 of block 0, just read that.
               // NOTE: This is prone to bugs, if you rearrange the storage of 
               // of the counters within a block, so beware.
    READ_RAW128 // Read the first 128 bits from the block.  See comment above.
  }; 

  enum struct FingerprintReadMethodEnum{
    READ_SIMPLE,
    READ_CROSS,
    READ_BYTE  // Special optimization for 8-bit fingerprints that are byte 
               // aligned
    // RAW reads don't make sense here. We don't statically know which atom 
    // that needs to be read.  It may make sense with 128-bit atoms, but 
    // my benchmarking showed 64-bit atoms to be faster.
  };

  enum struct FingerprintComparisonMethodEnum{
    VARIABLE_COUNT,
    FIXED_COUNT_AGGRESSIVE,
    SEMI_FIXED
  };

  enum struct ReductionMethodEnum{
    POP_CNT, // Must only use when counters fit into a single atom
    PARALLEL_REDUCE,
    NAIVE_FULL_EXCLUSIVE_SCAN,
  };

  enum struct OverflowTrackingArrayHashingMethodEnum{
    // Daniel Lemire's fast hashing method
    LEMIRE_FINGERPRINT_MULTIPLY,
    RAW_BUCKET_HASH,
    CLUSTERED_BUCKET_HASH,
  };

  enum struct InsertStatus{
    FAILED_TO_INSERT = 0,
    PLACED_IN_PRIMARY_BUCKET = 1,
    PLACED_IN_SECONDARY_BUCKET = 2
  };


} // End of CompressedCuckoo namespace

#endif
