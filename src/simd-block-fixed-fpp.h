// Copied from Apache Impala (incubating), usable under the terms in the Apache License,
// Version 2.0.

// This is a block Bloom filter (from Putze et al.'s "Cache-, Hash- and Space-Efficient
// Bloom Filters") with some twists:
//
// 1. Each block is a split Bloom filter - see Section 2.1 of Broder and Mitzenmacher's
// "Network Applications of Bloom Filters: A Survey".
//
// 2. The number of bits set per Add() is contant in order to take advantage of SIMD
// instructions.

#pragma once

#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <new>

#include <immintrin.h>

#include "hashutil.h"

using uint32_t = ::std::uint32_t;
using uint64_t = ::std::uint64_t;

__attribute__((always_inline))
inline uint32_t reduce(uint32_t hash, uint32_t n) {
    // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
    return (uint32_t) (((uint64_t) hash * n) >> 32);
}

inline uint64_t rotl64(uint64_t n, unsigned int c) {
    // assumes width is a power of 2
    const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
    // assert ( (c<=mask) &&"rotate by type width or more");
    c &= mask;
    return (n << c) | ( n >> ((-c) & mask));
}

template<typename HashFamily = ::cuckoofilter::TwoIndependentMultiplyShift>
class SimdBlockFilterFixed {
 private:
  // The filter is divided up into Buckets:
  using Bucket = uint32_t[8];

  const int bucketCount;

  Bucket* directory_;

  HashFamily hasher_;

 public:
  // Consumes at most (1 << log_heap_space) bytes on the heap:
  explicit SimdBlockFilterFixed(const int bits);
  ~SimdBlockFilterFixed() noexcept;
  void Add(const uint64_t key) noexcept;
  bool Find(const uint64_t key) const noexcept;
  uint64_t SizeInBytes() const { return sizeof(Bucket) * bucketCount; }

 private:
  // A helper function for Insert()/Find(). Turns a 32-bit hash into a 256-bit Bucket
  // with 1 single 1-bit set in each 32-bit lane.
  static __m256i MakeMask(const uint32_t hash) noexcept;

};

template<typename HashFamily>
SimdBlockFilterFixed<HashFamily>::SimdBlockFilterFixed(const int bits)
    // bits / 16: fpp 0.1777%, 75.1%
    // bits / 20: fpp 0.4384%, 63.4%
    // bits / 22: fpp 0.6692%, 61.1%
    // bits / 24: fpp 0.9765%, 59.7% <<== seems to be best (1% fpp seems important)
    // bits / 26: fpp 1.3769%, 59.3%
    // bits / 28: fpp 1.9197%, 60.3%
    // bits / 32: fpp 3.3280%, 63.0%
  : bucketCount(::std::max(1, bits / 24)),
    directory_(nullptr),
    hasher_() {
  if (!__builtin_cpu_supports("avx2")) {
    throw ::std::runtime_error("SimdBlockFilterFixed does not work without AVX2 instructions");
  }
  const size_t alloc_size = bucketCount * sizeof(Bucket);
  const int malloc_failed =
      posix_memalign(reinterpret_cast<void**>(&directory_), 64, alloc_size);
  if (malloc_failed) throw ::std::bad_alloc();
  memset(directory_, 0, alloc_size);
}

template<typename HashFamily>
SimdBlockFilterFixed<HashFamily>::~SimdBlockFilterFixed() noexcept {
  free(directory_);
  directory_ = nullptr;
}

// The SIMD reinterpret_casts technically violate C++'s strict aliasing rules. However, we
// compile with -fno-strict-aliasing.
template <typename HashFamily>
[[gnu::always_inline]] inline __m256i
SimdBlockFilterFixed<HashFamily>::MakeMask(const uint32_t hash) noexcept {
  const __m256i ones = _mm256_set1_epi32(1);
  // Odd contants for hashing:
  const __m256i rehash = _mm256_setr_epi32(0x47b6137bU, 0x44974d91U, 0x8824ad5bU,
      0xa2b7289dU, 0x705495c7U, 0x2df1424bU, 0x9efc4947U, 0x5c6bfb31U);
  // Load hash into a YMM register, repeated eight times
  __m256i hash_data = _mm256_set1_epi32(hash);
  // Multiply-shift hashing ala Dietzfelbinger et al.: multiply 'hash' by eight different
  // odd constants, then keep the 5 most significant bits from each product.
  hash_data = _mm256_mullo_epi32(rehash, hash_data);
  hash_data = _mm256_srli_epi32(hash_data, 27);
  // Use these 5 bits to shift a single bit to a location in each 32-bit lane
  return _mm256_sllv_epi32(ones, hash_data);
}

template <typename HashFamily>
[[gnu::always_inline]] inline void
SimdBlockFilterFixed<HashFamily>::Add(const uint64_t key) noexcept {
  const auto hash = hasher_(key);
  const uint32_t bucket_idx = reduce(rotl64(hash, 32), bucketCount);
  const __m256i mask = MakeMask(hash);
  __m256i* const bucket = &reinterpret_cast<__m256i*>(directory_)[bucket_idx];
  _mm256_store_si256(bucket, _mm256_or_si256(*bucket, mask));
}

template <typename HashFamily>
[[gnu::always_inline]] inline bool
SimdBlockFilterFixed<HashFamily>::Find(const uint64_t key) const noexcept {
  const auto hash = hasher_(key);
  const uint32_t bucket_idx = reduce(rotl64(hash, 32), bucketCount);
  const __m256i mask = MakeMask(hash);
  const __m256i bucket = reinterpret_cast<__m256i*>(directory_)[bucket_idx];
  // We should return true if 'bucket' has a one wherever 'mask' does. _mm256_testc_si256
  // takes the negation of its first argument and ands that with its second argument. In
  // our case, the result is zero everywhere iff there is a one in 'bucket' wherever
  // 'mask' is one. testc returns 1 if the result is 0 everywhere and returns 0 otherwise.
  return _mm256_testc_si256(bucket, mask);
}

