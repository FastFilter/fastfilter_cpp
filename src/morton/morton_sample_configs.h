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
// This file provides some sample configurations for Morton filters.  Given 
// that there are 16 template parameters, I wanted to hide this complexity
// from a standard user.

#ifndef _MORTON_SAMPLE_CONFIGS_H
#define _MORTON_SAMPLE_CONFIGS_H

#include "compressed_cuckoo_filter.h"
#include "fixed_point.h"

namespace CompressedCuckoo{

// Below are some sample configurations.  Note that the bits per
// item is different among the distinct configurations.  For example,
// Morton3_8 and Morton7_8 use about 11.7 and 9.98 bits per item, respectively,
// even though they use the same fingerprint length.  LF is short for load 
// factor in the table.

// Name       |  Description                           |  Bits per item (LF=0.95)
// Morton1_8; // 1-slot bucket with 8-bit fingerprints             11.7
// Morton3_8; // 3-slot buckets with 8-bit fingerprints            11.7
// Morton7_8; // 7-slot buckets with 8-bit fingerprints            9.98
// Morton15_8; // 15-slot buckets with 8-bit fingerprints          9.98


// Morton3_6; // 3-slot buckets with 6-bit fingerprints            8.84
// Morton7_6; // 7-slot buckets with 6-bit fingerprints            8.17
// Morton15_6; // 15-slot buckets with 6-bit fingerprints          8.84

// Morton3_12; // 3-slot buckets with 12-bit fingerprints          15.0
// Morton7_12; // 7-slot buckets with 12-bit fingerprints          13.5
// Morton15_12; // 15-slot buckets with 12-bit fingerprints        15.0

// Morton3_16; // 3-slot buckets with 16-bit fingerprints          20.0
// Morton7_16; // 7-slot buckets with 16-bit fingerprints          19.2
// Morton15_16; // 15-slot buckets with 16-bit fingerprints        18.6

// Morton3_18; // 3-slot buckets with 18-bit fingerprints          22.5
// Morton7_18; // 7-slot buckets with 18-bit fingerprints          21.6 
// Morton15_18; // 15-slot buckets with 18-bit fingerprints        20.7 

// TODO: Make some configurations with long fingerprints and rigorously test all

// Beginning of defintions

// Resizing is not enabled by default because it reduces performance, as 
// it requires additional computation that is not free.
// Set this to true to enable filters to self-resize
// Self-resizing impacts the false positive rate, as it reduces the effective 
// fingerprint length by one with each power-of-two doubling of the filter size.
// Be careful.
constexpr bool resizing_enabled = false; 

// Note: Using 1-slot buckets is much more brittle than using 3 slots or more due to an 
// increased likelihood of cycles disallowing new items to be inserted.  A stash 
// with tens to hundreds of items is necessary to make this work for filters with 
// hundreds of millions of elements that are filled to a FSA load ($\alpha_C$) of 0.95.

// 1-slot bucket with 8-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 11.7 bits per item at a load of 0.95
// Overflow Tracking Array: 16 x 1 bit
// Fullness Counters Array: 64 x 2 bits
// Fingerprint Storage Array: 46 x 8 bits
constexpr double target_compression_ratio_1_8 = 0.35937;
constexpr SerializedFixedPoint target_compression_ratio_sfp_1_8 = 
  FixedPoint(target_compression_ratio_1_8).serialize();
typedef CompressedCuckooFilter<
  1, // slots per bucket
  8, // fingerprint length in bits
  16, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_1_8, 
  CounterReadMethodEnum::READ_SIMPLE, 
  FingerprintReadMethodEnum::READ_SIMPLE,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::RAW_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton1_8;

// Main configuration from the VLDB'18 paper
// 3-slot bucket with 8-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 11.7 bits per item at a load factor of 0.95
// Overflow Tracking Array: 16 x 1 bit
// Fullness Counters Array: 64 x 2 bits
// Fingerprint Storage Array: 46 x 8 bits
constexpr double target_compression_ratio_3_8 = 0.23958;
constexpr SerializedFixedPoint target_compression_ratio_sfp_3_8 = 
  FixedPoint(target_compression_ratio_3_8).serialize();
typedef CompressedCuckooFilter<
  3, // slots per bucket
  8, // fingerprint length in bits
  16, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_3_8, 
  CounterReadMethodEnum::READ_SIMPLE, 
  FingerprintReadMethodEnum::READ_SIMPLE,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton3_8;

// 7-slot configuration from the VLDB'18 paper
// 7-slot bucket with 8-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 9.98 bits per item at a load factor of 0.95
// Overflow Tracking Array: 17 x 1 bit
// Fullness Counters Array: 21 x 3 bits
// Fingerprint Storage Array: 54 x 8 bits
constexpr double target_compression_ratio_7_8 = 0.36734;
constexpr SerializedFixedPoint target_compression_ratio_sfp_7_8 = 
  FixedPoint(target_compression_ratio_7_8).serialize();
typedef CompressedCuckooFilter<
  7, // slots per bucket
  8, // fingerprint length in bits
  17, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_7_8, 
  CounterReadMethodEnum::READ_RAW, 
  FingerprintReadMethodEnum::READ_SIMPLE,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton7_8;

// 15-slot configuration from the VLDB'18 paper
// 15-slot bucket with 8-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 9.98 bits per item at a load factor of 0.95
// Overflow Tracking Array: 16 x 1 bit
// Fullness Counters Array: 16 x 4 bits
// Fingerprint Storage Array: 54 x 8 bits
constexpr double target_compression_ratio_15_8 = 0.22499;
constexpr SerializedFixedPoint target_compression_ratio_sfp_15_8 = 
  FixedPoint(target_compression_ratio_15_8).serialize();
typedef CompressedCuckooFilter<
  15, // slots per bucket
  8, // fingerprint length in bits
  16, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_15_8, 
  CounterReadMethodEnum::READ_RAW, 
  FingerprintReadMethodEnum::READ_SIMPLE,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton15_8;

// 3-slot bucket with 6-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 8.84 bits per item at a load factor of 0.95
// Overflow Tracking Array: 16 x 1 bit
// Fullness Counters Array: 64 x 2 bits
// Fingerprint Storage Array: 61 x 6 bits
constexpr double target_compression_ratio_3_6 = 0.317708;
constexpr SerializedFixedPoint target_compression_ratio_sfp_3_6 = 
  FixedPoint(target_compression_ratio_3_6).serialize();
typedef CompressedCuckooFilter<
  3, // slots per bucket
  6, // fingerprint length in bits
  18, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_3_6, 
  CounterReadMethodEnum::READ_SIMPLE, 
  FingerprintReadMethodEnum::READ_CROSS,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton3_6;

// 7-slot bucket with 6-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 8.17 bits per item at a load factor of 0.95
// Overflow Tracking Array: 20 x 1 bit
// Fullness Counters Array: 32 x 3 bits
// Fingerprint Storage Array: 66 x 6 bits
constexpr double target_compression_ratio_7_6 = 0.294646;
constexpr SerializedFixedPoint target_compression_ratio_sfp_7_6 = 
  FixedPoint(target_compression_ratio_7_6).serialize();
typedef CompressedCuckooFilter<
  7, // slots per bucket
  6, // fingerprint length in bits
  20, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_7_6, 
  CounterReadMethodEnum::READ_CROSS, 
  FingerprintReadMethodEnum::READ_CROSS,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton7_6;

// FIXME: Too much memory use
// 15-slot bucket with 6-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 8.84 bits per item at a load factor of 0.95
// Overflow Tracking Array: 18 x 1 bit
// Fullness Counters Array: 32 x 4 bits
// Fingerprint Storage Array: 61 x 6 bits
constexpr double target_compression_ratio_15_6 = 0.127083;
constexpr SerializedFixedPoint target_compression_ratio_sfp_15_6 = 
  FixedPoint(target_compression_ratio_15_6).serialize();
typedef CompressedCuckooFilter<
  15, // slots per bucket
  6, // fingerprint length in bits
  18, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_15_6, 
  CounterReadMethodEnum::READ_SIMPLE, 
  FingerprintReadMethodEnum::READ_CROSS,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton15_6;

// 3-slot bucket with 12-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 15.0 bits per item at a load factor of 0.95
// Overflow Tracking Array: 16 x 1 bit
// Fullness Counters Array: 32 x 2 bits
// Fingerprint Storage Array: 36 x 12 bits
constexpr double target_compression_ratio_3_12 = 0.374999;
constexpr SerializedFixedPoint target_compression_ratio_sfp_3_12 = 
  FixedPoint(target_compression_ratio_3_12).serialize();
typedef CompressedCuckooFilter<
  3, // slots per bucket
  12, // fingerprint length in bits
  16, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_3_12, 
  CounterReadMethodEnum::READ_RAW, 
  FingerprintReadMethodEnum::READ_CROSS,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton3_12;

// 7-slot bucket with 12-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 13.5 bits per item at a load factor of 0.95
// Overflow Tracking Array: 16 x 1 bit
// Fullness Counters Array: 16 x 3 bits
// Fingerprint Storage Array: 40 x 12 bits
constexpr double target_compression_ratio_7_12 = 0.357142;
constexpr SerializedFixedPoint target_compression_ratio_sfp_7_12 = 
  FixedPoint(target_compression_ratio_7_12).serialize();
typedef CompressedCuckooFilter<
  7, // slots per bucket
  12, // fingerprint length in bits
  16, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_7_12, 
  CounterReadMethodEnum::READ_RAW, 
  FingerprintReadMethodEnum::READ_CROSS,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton7_12;

// 15-slot bucket with 12-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 15.0 bits per item at a load factor of 0.95
// Overflow Tracking Array: 16 x 1 bit
// Fullness Counters Array: 16 x 4 bits
// Fingerprint Storage Array: 36 x 12 bits
constexpr double target_compression_ratio_15_12 = 0.149999;
constexpr SerializedFixedPoint target_compression_ratio_sfp_15_12 = 
  FixedPoint(target_compression_ratio_15_12).serialize();
typedef CompressedCuckooFilter<
  15, // slots per bucket
  12, // fingerprint length in bits
  16, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_15_12, 
  CounterReadMethodEnum::READ_RAW, 
  FingerprintReadMethodEnum::READ_CROSS,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton15_12;

// 3-slot bucket with 8-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 20.0 bits per item at a load factor of 0.95
// Overflow Tracking Array: 16 x 1 bit
// Fullness Counters Array: 32 x 2 bits
// Fingerprint Storage Array: 27 x 16 bits
constexpr double target_compression_ratio_3_16 = 0.281249;
constexpr SerializedFixedPoint target_compression_ratio_sfp_3_16 = 
  FixedPoint(target_compression_ratio_3_16).serialize();
typedef CompressedCuckooFilter<
  3, // slots per bucket
  16, // fingerprint length in bits
  16, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_3_16, 
  CounterReadMethodEnum::READ_RAW, 
  FingerprintReadMethodEnum::READ_SIMPLE,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton3_16;

// 7-slot bucket with 16-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 19.24 bits per item at a load factor of 0.95
// Overflow Tracking Array: 16 x 1 bit
// Fullness Counters Array: 16 x 3 bits
// Fingerprint Storage Array: 28 x 16 bits
constexpr double target_compression_ratio_7_16 = 0.249999;
constexpr SerializedFixedPoint target_compression_ratio_sfp_7_16 = 
  FixedPoint(target_compression_ratio_7_16).serialize();
typedef CompressedCuckooFilter<
  7, // slots per bucket
  16, // fingerprint length in bits
  16, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_7_8, 
  CounterReadMethodEnum::READ_RAW, 
  FingerprintReadMethodEnum::READ_SIMPLE,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton7_16;

// 15-slot bucket with 16-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 18.6 bits per item at a load factor of 0.95
// Overflow Tracking Array: 16 x 1 bit
// Fullness Counters Array: 8 x 4 bits
// Fingerprint Storage Array: 29 x 16 bits
constexpr double target_compression_ratio_15_16 = 0.241666;
constexpr SerializedFixedPoint target_compression_ratio_sfp_15_16 = 
  FixedPoint(target_compression_ratio_15_16).serialize();
typedef CompressedCuckooFilter<
  15, // slots per bucket
  16, // fingerprint length in bits
  16, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_15_16, 
  CounterReadMethodEnum::READ_RAW, 
  FingerprintReadMethodEnum::READ_SIMPLE,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::LEMIRE_FINGERPRINT_MULTIPLY,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton15_16;

// 3-slot bucket with 18-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 22.5 bits per item at a load factor of 0.95
// Overflow Tracking Array: 16 x 1 bit
// Fullness Counters Array: 32 x 2 bits
// Fingerprint Storage Array: 24 x 18 bits
constexpr double target_compression_ratio_3_18 = 0.249999;
constexpr SerializedFixedPoint target_compression_ratio_sfp_3_18 =
  FixedPoint(target_compression_ratio_3_18).serialize();
typedef CompressedCuckooFilter<
  3, // slots per bucket
  18, // fingerprint length in bits
  16, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_3_18,
  CounterReadMethodEnum::READ_RAW,
  FingerprintReadMethodEnum::READ_CROSS,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET,
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton3_18;

// 7-slot bucket with 18-bit fingerprints
// parameter choices optimize for performance over minimizing storage costs.
// uses 21.6 bits per item at a load factor of 0.95
// overflow tracking array: 14 x 1 bit
// Fullness Counters Array: 16 x 3 bits
// Fingerprint Storage Array: 25 x 18 bits
constexpr double target_compression_ratio_7_18 = 0.223214;
constexpr SerializedFixedPoint target_compression_ratio_sfp_7_18 = 
  FixedPoint(target_compression_ratio_7_18).serialize();
typedef CompressedCuckooFilter<
  7, // slots per bucket
  18, // fingerprint length in bits
  14, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_7_18, 
  CounterReadMethodEnum::READ_RAW, 
  FingerprintReadMethodEnum::READ_CROSS,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton7_18;

// 15-slot bucket with 18-bit fingerprints
// Parameter choices optimize for performance over minimizing storage costs.
// Uses 20.7 bits per item at a load factor of 0.95
// Overflow Tracking Array: 12 x 1 bit
// Fullness Counters Array: 8 x 4 bits
// Fingerprint Storage Array: 26 x 18 bits
constexpr double target_compression_ratio_15_18 = 0.216666;
constexpr SerializedFixedPoint target_compression_ratio_sfp_15_18 = 
  FixedPoint(target_compression_ratio_15_18).serialize();
typedef CompressedCuckooFilter<
  15, // slots per bucket
  18, // fingerprint length in bits
  12, // overflow tracking array length in bits
  512, // block size in bits (should evenly divide into cache line)
  target_compression_ratio_sfp_15_18, 
  CounterReadMethodEnum::READ_RAW, 
  FingerprintReadMethodEnum::READ_SIMPLE,
  ReductionMethodEnum::POP_CNT,
  AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET, 
  OverflowTrackingArrayHashingMethodEnum::LEMIRE_FINGERPRINT_MULTIPLY,
  resizing_enabled, // resizing enabled
  true, // remapping of items from first bucket enabled
  true, // collision resolution enabled
  true, // Morton filter functionality enabled
  false, // Block fullness array enabled
  true,  // Handle conflicts on insertions enabled
  FingerprintComparisonMethodEnum::VARIABLE_COUNT
  > Morton15_18;



} // End of CompressedCuckoo namespace

#endif
