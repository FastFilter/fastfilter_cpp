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

template<typename HashFamily = ::cuckoofilter::TwoIndependentMultiplyShift>
class SimdBlockFilter {
 private:
  // The filter is divided up into Buckets:
  using Bucket = uint32_t[8];

  // log2(number of bytes in a bucket):
  static constexpr int LOG_BUCKET_BYTE_SIZE = 5;

  static_assert(
      (1 << LOG_BUCKET_BYTE_SIZE) == sizeof(Bucket) && sizeof(Bucket) == sizeof(__m256i),
      "Bucket sizing has gone awry.");

  // log_num_buckets_ is the log (base 2) of the number of buckets in the directory:
  const int log_num_buckets_;

  // directory_mask_ is (1 << log_num_buckets_) - 1. It is precomputed in the contructor
  // for efficiency reasons:
  const uint32_t directory_mask_;

  Bucket* directory_;

  HashFamily hasher_;

 public:
  // Consumes at most (1 << log_heap_space) bytes on the heap:
  explicit SimdBlockFilter(const int log_heap_space);
  SimdBlockFilter(SimdBlockFilter&& that)
    : log_num_buckets_(that.log_num_buckets_),
      directory_mask_(that.directory_mask_),
      directory_(that.directory_),
      hasher_(that.hasher_) {}
  ~SimdBlockFilter() noexcept;
  void Add(const uint64_t key) noexcept;
  bool Find(const uint64_t key) const noexcept;
  uint64_t SizeInBytes() const { return sizeof(Bucket) * (1ull << log_num_buckets_); }

 private:
  // A helper function for Insert()/Find(). Turns a 32-bit hash into a 256-bit Bucket
  // with 1 single 1-bit set in each 32-bit lane.
  static __m256i MakeMask(const uint32_t hash) noexcept;

  SimdBlockFilter(const SimdBlockFilter&) = delete;
  void operator=(const SimdBlockFilter&) = delete;
};

template<typename HashFamily>
SimdBlockFilter<HashFamily>::SimdBlockFilter(const int log_heap_space)
  :  // Since log_heap_space is in bytes, we need to convert it to the number of Buckets
     // we will use.
    log_num_buckets_(::std::max(1, log_heap_space - LOG_BUCKET_BYTE_SIZE)),
    // Don't use log_num_buckets_ if it will lead to undefined behavior by a shift that is
    // too large.
    directory_mask_((1ull << ::std::min(63, log_num_buckets_)) - 1),
    directory_(nullptr),
    hasher_() {
  if (!__builtin_cpu_supports("avx2")) {
    throw ::std::runtime_error("SimdBlockFilter does not work without AVX2 instructions");
  }
  const size_t alloc_size = 1ull << (log_num_buckets_ + LOG_BUCKET_BYTE_SIZE);
  const int malloc_failed =
      posix_memalign(reinterpret_cast<void**>(&directory_), 64, alloc_size);
  if (malloc_failed) throw ::std::bad_alloc();
  memset(directory_, 0, alloc_size);
}

template<typename HashFamily>
SimdBlockFilter<HashFamily>::~SimdBlockFilter() noexcept {
  free(directory_);
  directory_ = nullptr;
}

// The SIMD reinterpret_casts technically violate C++'s strict aliasing rules. However, we
// compile with -fno-strict-aliasing.
template <typename HashFamily>
[[gnu::always_inline]] inline __m256i
SimdBlockFilter<HashFamily>::MakeMask(const uint32_t hash) noexcept {
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
SimdBlockFilter<HashFamily>::Add(const uint64_t key) noexcept {
  const auto hash = hasher_(key);
  const uint32_t bucket_idx = hash & directory_mask_;
  const __m256i mask = MakeMask(hash >> log_num_buckets_);
  __m256i* const bucket = &reinterpret_cast<__m256i*>(directory_)[bucket_idx];
  _mm256_store_si256(bucket, _mm256_or_si256(*bucket, mask));
}

template <typename HashFamily>
[[gnu::always_inline]] inline bool
SimdBlockFilter<HashFamily>::Find(const uint64_t key) const noexcept {
  const auto hash = hasher_(key);
  const uint32_t bucket_idx = hash & directory_mask_;
  const __m256i mask = MakeMask(hash >> log_num_buckets_);
  const __m256i bucket = reinterpret_cast<__m256i*>(directory_)[bucket_idx];
  // We should return true if 'bucket' has a one wherever 'mask' does. _mm256_testc_si256
  // takes the negation of its first argument and ands that with its second argument. In
  // our case, the result is zero everywhere iff there is a one in 'bucket' wherever
  // 'mask' is one. testc returns 1 if the result is 0 everywhere and returns 0 otherwise.
  return _mm256_testc_si256(bucket, mask);
}




/// Rest is copied and pasted to work over 64-byte blocks

template <typename HashFamily = ::cuckoofilter::TwoIndependentMultiplyShift>
class SimdBlockFilter64 {
private:
  // The filter is divided up into Buckets:
  using Bucket = uint32_t[16];

  // log2(number of bytes in a bucket):
  static constexpr int LOG_BUCKET_BYTE_SIZE = 6;

  static_assert((1 << LOG_BUCKET_BYTE_SIZE) == sizeof(Bucket) &&
                    sizeof(Bucket) == 2 * sizeof(__m256i),
                "Bucket sizing has gone awry.");

  // log_num_buckets_ is the log (base 2) of the number of buckets in the
  // directory:
  const int log_num_buckets_;

  // directory_mask_ is (1 << log_num_buckets_) - 1. It is precomputed in the
  // contructor for efficiency reasons:
  const uint32_t directory_mask_;

  Bucket *directory_;

  HashFamily hasher_;

public:
  // Consumes at most (1 << log_heap_space) bytes on the heap:
  explicit SimdBlockFilter64(const int log_heap_space);
  SimdBlockFilter64(SimdBlockFilter64 &&that)
      : log_num_buckets_(that.log_num_buckets_),
        directory_mask_(that.directory_mask_), directory_(that.directory_),
        hasher_(that.hasher_) {}
  ~SimdBlockFilter64() noexcept;
  void Add(const uint64_t key) noexcept;
  bool Find(const uint64_t key) const noexcept;
  uint64_t SizeInBytes() const {
    return sizeof(Bucket) * (1ull << log_num_buckets_);
  }

private:
  static void MakeMask(const uint32_t hash, __m256i *out1,
                       __m256i *out2) noexcept;

  SimdBlockFilter64(const SimdBlockFilter64 &) = delete;
  void operator=(const SimdBlockFilter64 &) = delete;
};

template <typename HashFamily>
SimdBlockFilter64<HashFamily>::SimdBlockFilter64(const int log_heap_space)
    : // Since log_heap_space is in bytes, we need to convert it to the number
      // of Buckets we will use.
      log_num_buckets_(::std::max(1, log_heap_space - LOG_BUCKET_BYTE_SIZE)),
      // Don't use log_num_buckets_ if it will lead to undefined behavior by a
      // shift that is too large.
      directory_mask_((1ull << ::std::min(63, log_num_buckets_)) - 1),
      directory_(nullptr), hasher_() {
  if (!__builtin_cpu_supports("avx2")) {
    throw ::std::runtime_error(
        "SimdBlockFilter64 does not work without AVX2 instructions");
  }
  const size_t alloc_size = 1ull << (log_num_buckets_ + LOG_BUCKET_BYTE_SIZE);
  const int malloc_failed =
      posix_memalign(reinterpret_cast<void **>(&directory_), 64, alloc_size);
  if (malloc_failed)
    throw ::std::bad_alloc();
  memset(directory_, 0, alloc_size);
}

template <typename HashFamily>
SimdBlockFilter64<HashFamily>::~SimdBlockFilter64() noexcept {
  free(directory_);
  directory_ = nullptr;
}

// with AVX-512, this becomes a single instruction
static inline __m256i hacked_mm256_mullo_epi64(__m256i x, __m256i ml,
                                               __m256i mh) {
  __m256i xl = x;
  //      _mm256_and_si256(x, _mm256_set1_epi64x(UINT64_C(0x00000000ffffffff)));
  //  __m256i xh = _mm256_srli_epi64(x, 32);
  //__m256i hl = _mm256_slli_epi64(_mm256_mul_epu32(xh, ml), 32);
  __m256i lh = _mm256_slli_epi64(_mm256_mul_epu32(xl, mh), 32);
  __m256i ll = _mm256_mul_epu32(xl, ml);
  // return _mm256_add_epi64(ll, _mm256_add_epi64(hl, lh));
  return _mm256_add_epi64(lh, ll);
}

// The SIMD reinterpret_casts technically violate C++'s strict aliasing rules.
// However, we compile with -fno-strict-aliasing.
template <typename HashFamily>
[[gnu::always_inline]] inline void
SimdBlockFilter64<HashFamily>::MakeMask(const uint32_t hash, __m256i *out1,
                                        __m256i *out2) noexcept {
  const __m256i ones = _mm256_set1_epi64x(1);
  // Odd contants for hashing:
  const __m256i rehash1_l = _mm256_setr_epi64x(
      0x53214365047b6137 & 0xffffffff, 0x2c5635344974d91 & 0xffffffff,
      0x7fe299d78824ad5b & 0xffffffff, 0xc01ac48e4d29f115 & 0xffffffff);
  const __m256i rehash1_h = _mm256_setr_epi64x(
      UINT64_C(0x53214365047b6137) >> 32, UINT64_C(0x2c5635344974d91) >> 32,
      UINT64_C(0x7fe299d78824ad5b) >> 32, UINT64_C(0xc01ac48e4d29f115) >> 32);

  const __m256i rehash2_l = _mm256_setr_epi64x(
      0x7bdeb6734f95e2e3 & 0xffffffff, 0x2ec75a90a4e6ad3d & 0xffffffff,
      0x3d485cae00ae48fd & 0xffffffff, 0xe7d0f0c09b59d29b & 0xffffffff);
  const __m256i rehash2_h = _mm256_setr_epi64x(
      UINT64_C(0x7bdeb6734f95e2e3) >> 32, UINT64_C(0x2ec75a90a4e6ad3d) >> 32,
      UINT64_C(0x3d485cae00ae48fd) >> 32, UINT64_C(0xe7d0f0c09b59d29b) >> 32);

  __m256i hash_data = _mm256_set1_epi64x(hash);

  // Multiply-shift hashing ala Dietzfelbinger et al.: multiply 'hash' by eight
  // different odd constants, then keep the 6 most significant bits from each
  // product.
  __m256i hash_data1 = hacked_mm256_mullo_epi64(
      hash_data, rehash1_l,
      rehash1_h); //_mm256_mullo_epi64(rehash1, hash_data);
  __m256i hash_data2 = hacked_mm256_mullo_epi64(
      hash_data, rehash2_l,
      rehash2_h); //_mm256_mullo_epi64(rehash2, hash_data);

  hash_data1 = _mm256_and_si256(_mm256_srli_epi64(hash_data1, 32),
                                _mm256_set1_epi64x(63));
  hash_data2 = _mm256_and_si256(_mm256_srli_epi64(hash_data2, 32),
                                _mm256_set1_epi64x(63));
  // Use these 6 bits to shift a single bit to a location in each 32-bit lane
  *out1 = _mm256_sllv_epi64(ones, hash_data1);
  *out2 = _mm256_sllv_epi64(ones, hash_data2);
}


template <typename HashFamily>
[[gnu::always_inline]] inline void
SimdBlockFilter64<HashFamily>::Add(const uint64_t key) noexcept {
  const auto hash = hasher_(key);
  const uint32_t bucket_idx = hash & directory_mask_;
  __m256i mask1, mask2;
  MakeMask(hash >> log_num_buckets_, &mask1, &mask2);
  __m256i *const bucket1 =
      &reinterpret_cast<__m256i *>(directory_)[2 * bucket_idx];
  __m256i *const bucket2 =
      &reinterpret_cast<__m256i *>(directory_)[2 * bucket_idx + 1];

  _mm256_store_si256(bucket1, _mm256_or_si256(*bucket1, mask1));
  _mm256_store_si256(bucket2, _mm256_or_si256(*bucket2, mask2));
}

template <typename HashFamily>
[[gnu::always_inline]] inline bool
SimdBlockFilter64<HashFamily>::Find(const uint64_t key) const noexcept {
  const auto hash = hasher_(key);
  const uint32_t bucket_idx = hash & directory_mask_;
  __m256i mask1, mask2;
  MakeMask(hash >> log_num_buckets_, &mask1, &mask2);
  const __m256i bucket1 =
      reinterpret_cast<__m256i *>(directory_)[2 * bucket_idx];
  const __m256i bucket2 =
      reinterpret_cast<__m256i *>(directory_)[2 * bucket_idx + 1];
  return _mm256_testc_si256(bucket1, mask1) &
         _mm256_testc_si256(bucket2, mask2);
}
