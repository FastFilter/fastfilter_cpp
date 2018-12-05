#ifndef COUNTING_BLOOM_FILTER_BLOOM_FILTER_H_
#define COUNTING_BLOOM_FILTER_BLOOM_FILTER_H_

#include <algorithm>
#include <assert.h>

#include "hashutil.h"

#if defined(__BMI2__)
#include <immintrin.h>
#endif

using namespace std;
using namespace hashing;

namespace counting_bloomfilter {

inline int bitCount64(uint64_t x) {
    return __builtin_popcountll(x);
}

inline int select64(uint64_t x, int n) {
#if defined(__BMI2__)
    // with this, "add" is around 310 ns/key at 10000000 keys
    // from http://bitmagic.io/rank-select.html
    // https://github.com/Forceflow/libmorton/issues/6
    // This is a rather unusual usage of the pdep (bit deposit) operation,
    // as we use the x as the mask, and we use n as the value.
    // We deposit (move) the bits of 1 << n to the locations
    // defined by x.
    uint64_t d = _pdep_u64(1ULL << n, x);
    // and now we count the trailing zeroes, to find out
    // where the '1' was deposited
    return __builtin_ctzl(d);
    // return _tzcnt_u64(d);
#else
    // alternative implementation
    // with this, "add" is around 680 ns/key at 10000000 keys
    for(int i = 0; i < 64; i++) {
        if ((x & 1) == 1) {
            if (n-- == 0) {
                return i;
            }
        }
        x >>= 1;
    }
    return -1;
#endif
}

inline int numberOfLeadingZeros64(uint64_t x) {
    // If x is 0, the result is undefined.
    return __builtin_clzl(x);
}

enum Status {
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

inline uint32_t reduce(uint32_t hash, uint32_t n) {
  // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
  return (uint32_t)(((uint64_t)hash * n) >> 32);
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily = TwoIndependentMultiplyShift,
          int k = (int)((double)bits_per_item * 0.693147180559945 + 0.5)>
class CountingBloomFilter {

  uint64_t *data;
  size_t arrayLength;
  HashFamily hasher;

public:
  explicit CountingBloomFilter(const size_t n) : hasher() {
    size_t bitCount = 4 * n * bits_per_item;
    this->arrayLength = (bitCount + 63) / 64;
    data = new uint64_t[arrayLength]();
  }

  ~CountingBloomFilter() { delete[] data; }

  Status Add(const ItemType &item);

  Status Contain(const ItemType &item) const;

  size_t SizeInBytes() const { return arrayLength * 8; }
};

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
Status CountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    Add(const ItemType &key) {
  uint64_t hash = hasher(key);
  uint32_t a = (uint32_t)(hash >> 32);
  uint32_t b = (uint32_t)hash;
  for (int i = 0; i < k; i++) {
    uint index = reduce(a, this->arrayLength);
    data[index] += 1ULL << ((a << 2) & 63);
    a += b;
  }
  return Ok;
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
Status CountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    Contain(const ItemType &key) const {
  uint64_t hash = hasher(key);
  uint32_t a = (uint32_t)(hash >> 32);
  uint32_t b = (uint32_t)hash;
  for (int i = 0; i < k; i++) {
    uint index = reduce(a, this->arrayLength);
    if (((data[index] >> ((a << 2) & 63)) & 0xf) == 0) {
      return NotFound;
    }
    a += b;
  }
  return Ok;
}

// --------------------------------------------------------------------------------------

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily = TwoIndependentMultiplyShift,
          int k = (int)((double)bits_per_item * 0.693147180559945 + 0.5)>
class SuccinctCountingBloomFilter {

  uint64_t *data;
  uint64_t *counts;
  uint64_t *overflow;
  // uint8_t *realCount;
  size_t arrayLength;
  size_t overflowLength;
  size_t nextFreeOverflow;
  HashFamily hasher;

  void Increment(size_t group, int bit);
  int ReadCount(size_t group, int bit);

public:
  explicit SuccinctCountingBloomFilter(const size_t n) : hasher() {
    size_t bitCount = n * bits_per_item;
    this->arrayLength = (bitCount + 63) / 64;
    this->overflowLength = 100 + arrayLength / 100 * 12;
    data = new uint64_t[arrayLength]();
    counts = new uint64_t[arrayLength]();
    overflow = new uint64_t[overflowLength]();
    // realCount = new uint8_t[arrayLength * 64]();
    nextFreeOverflow = 0;
    for (size_t i = 0; i < overflowLength; i += 4) {
        overflow[i] = i + 4;
    }
  }

  ~SuccinctCountingBloomFilter() { delete[] data; delete[] counts; delete[] overflow; }

  Status Add(const ItemType &item);

  Status Contain(const ItemType &item) const;

  size_t SizeInBytes() const { return arrayLength * 8 * 2 + overflowLength * 8; }
};

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
Status SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    Add(const ItemType &key) {
  uint64_t hash = hasher(key);
  uint32_t a = (uint32_t)(hash >> 32);
  uint32_t b = (uint32_t)hash;
  for (int i = 0; i < k; i++) {
    uint group = reduce(a, this->arrayLength);
    Increment(group, a & 63);
    a += b;
  }
  return Ok;
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
void SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    Increment(size_t group, int bit) {
    // realCount[(group << 6) + bit]++;
    uint64_t m = data[group];
    uint64_t c = counts[group];
    if ((c & 0xc000000000000000L) != 0) {
        // an overflow entry, or overflowing now
        size_t index;
        if ((c & 0x8000000000000000L) == 0) {
            // convert to an overflow entry
            // allocate overflow
            index = nextFreeOverflow;
            if (index >= overflowLength) {
                ::std::cout << "ERROR: overflow too small\n";
                data[group] |= 1ULL << bit;
                return;
            }
            nextFreeOverflow = (size_t) overflow[index];
            overflow[index] = 0;
            overflow[index + 1] = 0;
            overflow[index + 2] = 0;
            overflow[index + 3] = 0;
            // convert to a pointer
            for (int i = 0; i < 64; i++) {
                int n = ReadCount(group, i);
                overflow[index + i / 16] += n * (1ULL << (i * 4));
            }
            uint64_t count = 64;
            c = 0x8000000000000000L | (count << 32) | index;
            counts[group] = c;
        } else {
            // already
            index = (size_t) (c & 0x0fffffff);
            c += 1ULL << 32;
            counts[group] = c;
        }
        overflow[index + bit / 16] += (1ULL << (bit * 4));
        data[group] |= 1ULL << bit;
        return;
    }
    data[group] |= 1ULL << bit;
    int bitsBefore = bitCount64(m & (0xffffffffffffffffL >> (63 - bit)));
    int before = select64((c << 1) | 1, bitsBefore);
    int d = (m >> bit) & 1;
    int insertAt = before - d;
    uint64_t mask = (1ULL << insertAt) - 1;
    uint64_t left = c & ~mask;
    uint64_t right = c & mask;
    c = (left << 1) | ((1ULL ^ d) << insertAt) | right;
    counts[group] = c;
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
int SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    ReadCount(size_t group, int bit) {
    uint64_t m = data[group];
    uint64_t d = (m >> bit) & 1;
    if (d == 0) {
        return 0;
    }
    uint64_t c = counts[group];
    if ((c & 0x8000000000000000L) != 0) {
        size_t index = (size_t) (c & 0x0fffffff);
        uint64_t n = overflow[index + bit / 16];
        n >>= 4 * (bit & 0xf);
        return (int) (n & 15);
    }
    int bitsBefore = bitCount64(m & (0xffffffffffffffffL >> (63 - bit)));
    int bitPos = select64(c, bitsBefore - 1);
    uint64_t y = ((c << (63 - bitPos)) << 1) | (1ULL << (63 - bitPos));
    return numberOfLeadingZeros64(y) + 1;
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
Status SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    Contain(const ItemType &key) const {
  uint64_t hash = hasher(key);
  uint32_t a = (uint32_t)(hash >> 32);
  uint32_t b = (uint32_t)hash;
  for (int i = 0; i < k; i++) {
    uint group = reduce(a, this->arrayLength);
    if (((data[group] >> (a & 63)) & 1) == 0) {
      return NotFound;
    }
    a += b;
  }
  return Ok;
}

}
#endif