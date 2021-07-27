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

#define bitCount64(x) __builtin_popcountll(x)

#if defined(__BMI2__)

// use a macro, to ensure it is inlined
#define select64(A, B) _tzcnt_u64(_pdep_u64(1ULL << (B), (A)))

#else

/*
inline int select64(uint64_t x, int n) {
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
}
*/

/*

From https://github.com/splatlab/Sux
https://github.com/splatlab/Sux/blob/master/testselect64.cpp

*/


#define ONES_STEP_4 0x1111111111111111ULL
#define ONES_STEP_8 0x0101010101010101ULL
#define MSBS_STEP_8 (0x80L * ONES_STEP_8)

static int8_t SELECT_IN_BYTE[] = {
    -1, 0, 1, 0, 2, 0, 1, 0, 3,
    0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1,
    0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2,
    0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1,
    0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1,
    0, 3, 0, 1, 0, 2, 0, 1, 0, 7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2,
    0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1,
    0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3,
    0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1,
    0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2,
    0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1,
    0, 2, 0, 1, 0, -1, -1, -1, 1, -1, 2, 2, 1, -1, 3, 3, 1, 3, 2, 2, 1,
    -1, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, -1, 5, 5, 1, 5, 2,
    2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1,
    3, 2, 2, 1, -1, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4,
    4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1,
    5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2,
    2, 1, -1, 7, 7, 1, 7, 2, 2, 1, 7, 3, 3, 1, 3, 2, 2, 1, 7, 4, 4, 1,
    4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 7, 5, 5, 1, 5, 2, 2, 1, 5, 3,
    3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    7, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1, 4, 2,
    2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1,
    3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, -1, -1,
    -1, -1, -1, -1, -1, 2, -1, -1, -1, 3, -1, 3, 3, 2, -1, -1, -1, 4,
    -1, 4, 4, 2, -1, 4, 4, 3, 4, 3, 3, 2, -1, -1, -1, 5, -1, 5, 5, 2,
    -1, 5, 5, 3, 5, 3, 3, 2, -1, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3,
    3, 2, -1, -1, -1, 6, -1, 6, 6, 2, -1, 6, 6, 3, 6, 3, 3, 2, -1, 6,
    6, 4, 6, 4, 4, 2, 6, 4, 4, 3, 4, 3, 3, 2, -1, 6, 6, 5, 6, 5, 5, 2,
    6, 5, 5, 3, 5, 3, 3, 2, 6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3,
    3, 2, -1, -1, -1, 7, -1, 7, 7, 2, -1, 7, 7, 3, 7, 3, 3, 2, -1, 7,
    7, 4, 7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2, -1, 7, 7, 5, 7, 5, 5, 2,
    7, 5, 5, 3, 5, 3, 3, 2, 7, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3,
    3, 2, -1, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2, 7, 6, 6, 4,
    6, 4, 4, 2, 6, 4, 4, 3, 4, 3, 3, 2, 7, 6, 6, 5, 6, 5, 5, 2, 6, 5,
    5, 3, 5, 3, 3, 2, 6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1,
    -1, -1, -1, -1, -1, -1, 4, -1, -1, -1, 4, -1, 4, 4, 3, -1, -1, -1,
    -1, -1, -1, -1, 5, -1, -1, -1, 5, -1, 5, 5, 3, -1, -1, -1, 5, -1,
    5, 5, 4, -1, 5, 5, 4, 5, 4, 4, 3, -1, -1, -1, -1, -1, -1, -1, 6,
    -1, -1, -1, 6, -1, 6, 6, 3, -1, -1, -1, 6, -1, 6, 6, 4, -1, 6, 6,
    4, 6, 4, 4, 3, -1, -1, -1, 6, -1, 6, 6, 5, -1, 6, 6, 5, 6, 5, 5, 3,
    -1, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, -1, -1, -1, -1,
    -1, -1, -1, 7, -1, -1, -1, 7, -1, 7, 7, 3, -1, -1, -1, 7, -1, 7, 7,
    4, -1, 7, 7, 4, 7, 4, 4, 3, -1, -1, -1, 7, -1, 7, 7, 5, -1, 7, 7,
    5, 7, 5, 5, 3, -1, 7, 7, 5, 7, 5, 5, 4, 7, 5, 5, 4, 5, 4, 4, 3, -1,
    -1, -1, 7, -1, 7, 7, 6, -1, 7, 7, 6, 7, 6, 6, 3, -1, 7, 7, 6, 7, 6,
    6, 4, 7, 6, 6, 4, 6, 4, 4, 3, -1, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5,
    6, 5, 5, 3, 7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 5, -1, -1, -1, -1, -1,
    -1, -1, 5, -1, -1, -1, 5, -1, 5, 5, 4, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, 6, -1, -1, -1, -1, -1, -1, -1, 6,
    -1, -1, -1, 6, -1, 6, 6, 4, -1, -1, -1, -1, -1, -1, -1, 6, -1, -1,
    -1, 6, -1, 6, 6, 5, -1, -1, -1, 6, -1, 6, 6, 5, -1, 6, 6, 5, 6, 5,
    5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    7, -1, -1, -1, -1, -1, -1, -1, 7, -1, -1, -1, 7, -1, 7, 7, 4, -1,
    -1, -1, -1, -1, -1, -1, 7, -1, -1, -1, 7, -1, 7, 7, 5, -1, -1, -1,
    7, -1, 7, 7, 5, -1, 7, 7, 5, 7, 5, 5, 4, -1, -1, -1, -1, -1, -1,
    -1, 7, -1, -1, -1, 7, -1, 7, 7, 6, -1, -1, -1, 7, -1, 7, 7, 6, -1,
    7, 7, 6, 7, 6, 6, 4, -1, -1, -1, 7, -1, 7, 7, 6, -1, 7, 7, 6, 7, 6,
    6, 5, -1, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 4, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, 5, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, 6, -1, -1, -1, -1, -1, -1, -1, 6, -1, -1,
    -1, 6, -1, 6, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 7, -1, -1, -1, -1, -1, -1, -1, 7, -1, -1, -1, 7, -1, 7, 7, 5,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 7, -1,
    -1, -1, -1, -1, -1, -1, 7, -1, -1, -1, 7, -1, 7, 7, 6, -1, -1, -1,
    -1, -1, -1, -1, 7, -1, -1, -1, 7, -1, 7, 7, 6, -1, -1, -1, 7, -1,
    7, 7, 6, -1, 7, 7, 6, 7, 6, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 7, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 7, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 7, -1, -1, -1, -1, -1, -1,
    -1, 7, -1, -1, -1, 7, -1, 7, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, 7
};

inline int select64(uint64_t x, int n) {
    uint64_t byteSums = x - ((x & 0xa * ONES_STEP_4) >> 1);
    byteSums = (byteSums & 3 * ONES_STEP_4) +
            ((byteSums >> 2) & 3 * ONES_STEP_4);
    byteSums = (byteSums + (byteSums >> 4)) & 0x0f * ONES_STEP_8;
    byteSums *= ONES_STEP_8;
    // Phase 2: compare each byte sum with rank to obtain the relevant byte
    uint64_t rankStep8 = n * ONES_STEP_8;
    int byteOffset = (int) (((((rankStep8 | MSBS_STEP_8) - byteSums) & MSBS_STEP_8) >> 7) *
            ONES_STEP_8 >> 53) &
            ~0x7;
    int byteRank = (int) (n - (((byteSums << 8) >> byteOffset) & 0xFF));
    return byteOffset +
            SELECT_IN_BYTE[(int) ((x >> byteOffset) & 0xFF) | byteRank << 8];
}

#endif

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

// CountingBloomFilter --------------------------------------------------------------------------------------

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily = SimpleMixSplit,
          int k = (int)((double)bits_per_item * 0.693147180559945 + 0.5)>
class CountingBloomFilter {

  uint64_t *data;
  size_t arrayLength;
  HashFamily hasher;
  const int blockShift = 14;
  const int blockLen = 1 << blockShift;

  void AddBlock(uint32_t *tmp, int block, int len);

public:
  explicit CountingBloomFilter(const size_t n) : hasher() {
    size_t bitCount = 4 * n * bits_per_item;
    this->arrayLength = (bitCount + 63) / 64;
    data = new uint64_t[arrayLength]();
  }
  ~CountingBloomFilter() { delete[] data; }
  Status Add(const ItemType &item);
  Status AddAll(const vector<ItemType>& data, const size_t start, const size_t end);
  Status Remove(const ItemType &item);
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
    data[index] += 1ULL << ((a << 2) & 0x3f);
    a += b;
  }
  return Ok;
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
void CountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    AddBlock(uint32_t *tmp, int block, int len) {
  for (int i = 0; i < len; i++) {
    int index = tmp[(block << blockShift) + i];
    data[index >> 4] += 1ULL << ((index << 2) & 0x3f);
  }
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
Status CountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    AddAll(const vector<ItemType>& keys, const size_t start, const size_t end) {
  int blocks = 1 + arrayLength / blockLen;
  uint32_t *tmp = new uint32_t[blocks * blockLen];
  int *tmpLen = new int[blocks]();
  for (size_t i = start; i < end; i++) {
    uint64_t key = keys[i];
    uint64_t hash = hasher(key);
    uint32_t a = (uint32_t)(hash >> 32);
    uint32_t b = (uint32_t)hash;
    for (int j = 0; j < k; j++) {
      int index = reduce(a, this->arrayLength);
      int block = index >> blockShift;
      int len = tmpLen[block];
      tmp[(block << blockShift) + len] = (index << 4) + (a & 0xf);
      tmpLen[block] = len + 1;
      if (len + 1 == blockLen) {
        AddBlock(tmp, block, len + 1);
        tmpLen[block] = 0;
      }
      a += b;
    }
  }
  for (int block = 0; block < blocks; block++) {
    AddBlock(tmp, block, tmpLen[block]);
  }
  delete[] tmp;
  delete[] tmpLen;
  return Ok;
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
Status CountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    Remove(const ItemType &key) {
  uint64_t hash = hasher(key);
  uint32_t a = (uint32_t)(hash >> 32);
  uint32_t b = (uint32_t)hash;
  for (int i = 0; i < k; i++) {
    uint index = reduce(a, this->arrayLength);
    data[index] -= 1ULL << ((a << 2) & 0x3f);
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
    if (((data[index] >> ((a << 2) & 0x3f)) & 0xf) == 0) {
      return NotFound;
    }
    a += b;
  }
  return Ok;
}

// SuccinctCountingBloomFilter --------------------------------------------------------------------------------------

// #define VERIFY_COUNT

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily = SimpleMixSplit,
          int k = (int)((double)bits_per_item * 0.693147180559945 + 0.5)>
class SuccinctCountingBloomFilter {

  uint64_t *data;
  uint64_t *counts;
  uint64_t *overflow;
#ifdef VERIFY_COUNT
  uint8_t *realCount;
#endif
  size_t arrayLength;
  size_t overflowLength;
  size_t nextFreeOverflow;
  HashFamily hasher;
  const int blockShift = 14;
  const int blockLen = 1 << blockShift;

  void Increment(size_t group, int bit);
  void Decrement(size_t group, int bit);
  int ReadCount(size_t group, int bit);
  void AddBlock(uint32_t *tmp, int block, int len);

public:
  explicit SuccinctCountingBloomFilter(const size_t n) : hasher() {
    size_t bitCount = n * bits_per_item;
    this->arrayLength = (bitCount + 63) / 64;
    this->overflowLength = 100 + arrayLength / 100 * 12;
    data = new uint64_t[arrayLength]();
    counts = new uint64_t[arrayLength]();
    overflow = new uint64_t[overflowLength]();
#ifdef VERIFY_COUNT
    realCount = new uint8_t[arrayLength * 64]();
#endif
    nextFreeOverflow = 0;
    for (size_t i = 0; i < overflowLength; i += 4) {
        overflow[i] = i + 4;
    }
  }
  ~SuccinctCountingBloomFilter() { delete[] data; delete[] counts; delete[] overflow; }
  Status Add(const ItemType &item);
  Status AddAll(const vector<ItemType>& data, const size_t start, const size_t end);
  Status Remove(const ItemType &item);
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
    AddBlock(uint32_t *tmp, int block, int len) {
  for (int i = 0; i < len; i++) {
    uint32_t index = tmp[(block << blockShift) + i];
    uint32_t group = index >> 6;
    Increment(group, index & 63);
  }
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
Status SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    Remove(const ItemType &key) {
  uint64_t hash = hasher(key);
  uint32_t a = (uint32_t)(hash >> 32);
  uint32_t b = (uint32_t)hash;
  for (int i = 0; i < k; i++) {
    uint group = reduce(a, this->arrayLength);
    Decrement(group, a & 63);
    a += b;
  }
  return Ok;
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
Status SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    AddAll(const vector<ItemType>& keys, const size_t start, const size_t end) {
  int blocks = 1 + arrayLength / blockLen;
  uint32_t *tmp = new uint32_t[blocks * blockLen];
  int *tmpLen = new int[blocks]();
  for (size_t i = start; i < end; i++) {
    uint64_t key = keys[i];
    uint64_t hash = hasher(key);
    uint32_t a = (uint32_t)(hash >> 32);
    uint32_t b = (uint32_t)hash;
    for (int j = 0; j < k; j++) {
      int index = reduce(a, this->arrayLength);
      int block = index >> blockShift;
      int len = tmpLen[block];
      tmp[(block << blockShift) + len] = (index << 6) + (a & 63);
      tmpLen[block] = len + 1;
      if (len + 1 == blockLen) {
        AddBlock(tmp, block, len + 1);
        tmpLen[block] = 0;
      }
      a += b;
    }
  }
  for (int block = 0; block < blocks; block++) {
    AddBlock(tmp, block, tmpLen[block]);
  }
  delete[] tmp;
  delete[] tmpLen;
  return Ok;
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
void SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    Increment(size_t group, int bit) {
#ifdef VERIFY_COUNT
    realCount[(group << 6) + bit]++;
#endif
    uint64_t m = data[group];
    uint64_t c = counts[group];
    if ((c & 0xc000000000000000ULL) != 0) {
        // an overflow entry, or overflowing now
        size_t index;
        if ((c & 0x8000000000000000ULL) == 0) {
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
            c = 0x8000000000000000ULL | (count << 32) | index;
            counts[group] = c;
        } else {
            // already
            index = (size_t) (c & 0x0fffffffULL);
            c += 1ULL << 32;
            counts[group] = c;
        }
        overflow[index + bit / 16] += (1ULL << (bit * 4));
        data[group] |= 1ULL << bit;
    } else {
        data[group] |= 1ULL << bit;
        int bitsBefore = bitCount64(m & (0xffffffffffffffffULL >> (63 - bit)));
        int before = select64((c << 1) | 1, bitsBefore);
        int d = (m >> bit) & 1;
        int insertAt = before - d;
        uint64_t mask = (1ULL << insertAt) - 1;
        uint64_t left = c & ~mask;
        uint64_t right = c & mask;
        c = (left << 1) | ((1ULL ^ d) << insertAt) | right;
        counts[group] = c;
    }
#ifdef VERIFY_COUNT
    for(int b = 0; b < 64; b++) {
        if (realCount[(group << 6) + b] != ReadCount(group, b)) {
            ::std::cout << "group " << group << "/" << b << " of " << bit << "\n";
        }
    }
#endif
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
    if ((c & 0x8000000000000000ULL) != 0) {
        size_t index = (size_t) (c & 0x0fffffffULL);
        uint64_t n = overflow[index + bit / 16];
        n >>= 4 * (bit & 0xf);
        return (int) (n & 15);
    }
    int bitsBefore = bitCount64(m & (0xffffffffffffffffULL >> (63 - bit)));
    int bitPos = select64(c, bitsBefore - 1);
    uint64_t y = ((c << (63 - bitPos)) << 1) | (1ULL << (63 - bitPos));
    return numberOfLeadingZeros64(y) + 1;
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
void SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::
    Decrement(size_t group, int bit) {
#ifdef VERIFY_COUNT
    realCount[(group << 6) + bit]--;
#endif
    uint64_t m = data[group];
    uint64_t c = counts[group];
    if ((c & 0x8000000000000000ULL) != 0) {
        // an overflow entry
        size_t index = (size_t) (c & 0x0fffffffULL);
        size_t count = (size_t) (c >> 32) & 0x0fffffffULL;
        c -= 1ULL << 32;
        counts[group] = c;
        uint64_t n = overflow[index + bit / 16];
        overflow[index + bit / 16] = n - (1ULL << (bit * 4));
        n >>= 4 * (bit & 0xf);
        if ((n & 0xf) == 1) {
            data[group] &= ~(1ULL << bit);
        }
        if (count < 64) {
            // convert back to an inline entry, and free up the overflow entry
            uint64_t c2 = 0;
            for (int j = 63; j >= 0; j--) {
                int cj = (int) ((overflow[index + j / 16] >> (4 * j)) & 0xf);
                if (cj > 0) {
                    c2 = ((c2 << 1) | 1) << (cj - 1);
                }
            }
            counts[group] = c2;
            // free overflow
            overflow[index] = nextFreeOverflow;
            nextFreeOverflow = index;
        }
    } else {
        int bitsBefore = bitCount64(m & (0xffffffffffffffffULL >> (63 - bit)));
        int before = select64((c << 1) | 1, bitsBefore) - 1;
        int removeAt = max(0, before - 1);
        // remove the bit from the counter
        uint64_t mask = (1ULL << removeAt) - 1;
        uint64_t left = (c >> 1) & ~mask;
        uint64_t right= c & mask;
        counts[group] = left | right;
        uint64_t removed = (c >> removeAt) & 1;
        // possibly reset the data bit
        data[group] = m & ~(removed << bit);
    }
#ifdef VERIFY_COUNT
    for(int b = 0; b < 64; b++) {
        if (realCount[(group << 6) + b] != ReadCount(group, b)) {
            ::std::cout << "group- " << group << "/" << b << " of " << bit << "\n";
        }
    }
#endif
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

// SuccinctCountingBlockedBloomFilter --------------------------------------------------------------------------------------

// #define VERIFY_COUNT

template <typename ItemType, size_t bits_per_item, typename HashFamily,
          int k = (int)((double)bits_per_item * 0.693147180559945 + 0.5)>
class SuccinctCountingBlockedBloomFilter {
private:
  const int bucketCount;
  HashFamily hasher;
  uint64_t *data;
  uint64_t *counts;
  uint64_t *overflow;
  size_t overflowLength;
  size_t nextFreeOverflow;
#ifdef VERIFY_COUNT
  uint8_t *realCount;
#endif

  void Increment(size_t group, int bit);
  void Decrement(size_t group, int bit);
  int ReadCount(size_t group, int bit);
#ifdef VERIFY_COUNT
  void VerifyCount(size_t group, int bit, int line);
#endif

public:
  explicit SuccinctCountingBlockedBloomFilter(const int capacity);
  ~SuccinctCountingBlockedBloomFilter() noexcept;
  void Add(const uint64_t key) noexcept;
  void Remove(const uint64_t key) noexcept;
  bool Contain(const uint64_t key) const noexcept;
  uint64_t SizeInBytes() const {
      return 2 * 64 * bucketCount + 8 * overflowLength;
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily, k>::
    SuccinctCountingBlockedBloomFilter(const int capacity)
    : bucketCount(capacity * bits_per_item / 512), hasher() {
  const size_t alloc_size = bucketCount * (512 / 8);
  const int malloc_failed =
      posix_memalign(reinterpret_cast<void **>(&data), 64, alloc_size);
  if (malloc_failed)
    throw ::std::bad_alloc();
  memset(data, 0, alloc_size);
  size_t arrayLength = bucketCount * 8;
  overflowLength = 100 + arrayLength / 100 * 36;
  counts = new uint64_t[arrayLength]();
  overflow = new uint64_t[overflowLength]();
#ifdef VERIFY_COUNT
  realCount = new uint8_t[arrayLength * 64]();
#endif
  nextFreeOverflow = 0;
  for (size_t i = 0; i < overflowLength; i += 8) {
      overflow[i] = i + 8;
  }
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily, k>::
    ~SuccinctCountingBlockedBloomFilter() noexcept {
  free(data);
  delete[] counts;
  delete[] overflow;
}

static inline uint64_t rotl64(uint64_t n, unsigned int c) {
  // assumes width is a power of 2
  const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
  // assert ( (c<=mask) &&"rotate by type width or more");
  c &= mask;
  return (n << c) | (n >> ((-c) & mask));
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
void SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily, k>::
    Add(const uint64_t key) noexcept {
  const auto hash = hasher(key);
  const uint32_t bucket_start = reduce(rotl64(hash, 32), bucketCount) * 8;
  uint32_t a = (uint32_t)hash;
  if (k >= 3) {
    Increment(bucket_start + ((a >> 0) & 7), (a >> 3) & 0x3f);
    Increment(bucket_start + ((a >> 9) & 7), (a >> 12) & 0x3f);
    Increment(bucket_start + ((a >> 18) & 7), (a >> 21) & 0x3f);
//    data[bucket_start + ((a >> 0) & 7)] |= 1ULL << ((a >> 3) & 0x3f);
//    data[bucket_start + ((a >> 9) & 7)] |= 1ULL << ((a >> 12) & 0x3f);
//    data[bucket_start + ((a >> 18) & 7)] |= 1ULL << ((a >> 21) & 0x3f);
  }
  uint32_t b = (uint32_t)(hash >> 32);
  for (int i = 3; i < k; i++) {
    a += b;
    Increment(bucket_start + (a & 7), (a >> 3) & 0x3f);
//    data[bucket_start + (a & 7)] |= 1ULL << ((a >> 3) & 0x3f);
  }
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
void SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily, k>::
    Increment(size_t group, int bit) {
#ifdef VERIFY_COUNT
    realCount[(group << 6) + bit]++;
#endif
    uint64_t m = data[group];
    uint64_t c = counts[group];
    if ((c & 0xc000000000000000ULL) != 0) {
        // an overflow entry, or overflowing now
        size_t index;
        if ((c & 0x8000000000000000ULL) == 0) {
            // convert to an overflow entry
            // allocate overflow
            index = nextFreeOverflow;
            if (index >= overflowLength) {
                ::std::cout << "ERROR: overflow too small\n";
                data[group] |= 1ULL << bit;
                return;
            }
            nextFreeOverflow = (size_t) overflow[index];
            for (int i = 0; i < 8; i++) {
                overflow[index + i] = 0;
            }
            // convert to a pointer
            for (int i = 0; i < 64; i++) {
                int n = ReadCount(group, i);
                overflow[index + i / 8] += n * (1ULL << (i * 8));
            }
            uint64_t count = 64;
            c = 0x8000000000000000ULL | (count << 32) | index;
            counts[group] = c;
        } else {
            // already
            index = (size_t) (c & 0x0fffffffULL);
            c += 1ULL << 32;
            counts[group] = c;
        }
        overflow[index + bit / 8] += (1ULL << (bit * 8));
        data[group] |= 1ULL << bit;
    } else {
        data[group] |= 1ULL << bit;
        int bitsBefore = bitCount64(m & (0xffffffffffffffffULL >> (63 - bit)));
        int before = select64((c << 1) | 1, bitsBefore);
        int d = (m >> bit) & 1;
        int insertAt = before - d;
        uint64_t mask = (1ULL << insertAt) - 1;
        uint64_t left = c & ~mask;
        uint64_t right = c & mask;
        c = (left << 1) | ((1ULL ^ d) << insertAt) | right;
        counts[group] = c;
    }
#ifdef VERIFY_COUNT
    VerifyCount(group, bit, __LINE__);
#endif
}

#ifdef VERIFY_COUNT
template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
void SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily, k>::
    VerifyCount(size_t group, int bit, int line) {
    for(int b = 0; b < 64; b++) {
        if (realCount[(group << 6) + b] != ReadCount(group, b)) {
            ::std::cout << "group " << group << " bit " << b <<
            " expected " << (int) realCount[(group << 6) + b] <<
            " got " << ReadCount(group, b) <<
            " at bit " << bit << " at line " << line << "\n";
        }
    }
}
#endif

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
int SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily, k>::
    ReadCount(size_t group, int bit) {
    uint64_t m = data[group];
    uint64_t d = (m >> bit) & 1;
    if (d == 0) {
        return 0;
    }
    uint64_t c = counts[group];
    if ((c & 0x8000000000000000ULL) != 0) {
        size_t index = (size_t) (c & 0x0fffffffULL);
        uint64_t n = overflow[index + bit / 8];
        n >>= 8 * (bit & 0xff);
        return (int) (n & 0xff);
    }
    int bitsBefore = bitCount64(m & (0xffffffffffffffffULL >> (63 - bit)));
    int bitPos = select64(c, bitsBefore - 1);
    uint64_t y = ((c << (63 - bitPos)) << 1) | (1ULL << (63 - bitPos));
    return numberOfLeadingZeros64(y) + 1;
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
void SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily, k>::
    Remove(const uint64_t key) noexcept {
  const auto hash = hasher(key);
  const uint32_t bucket_start = reduce(rotl64(hash, 32), bucketCount) * 8;
  uint32_t a = (uint32_t)hash;
  if (k >= 3) {
    Decrement(bucket_start + ((a >> 0) & 7), (a >> 3) & 0x3f);
    Decrement(bucket_start + ((a >> 9) & 7), (a >> 12) & 0x3f);
    Decrement(bucket_start + ((a >> 18) & 7), (a >> 21) & 0x3f);
  }
  uint32_t b = (uint32_t)(hash >> 32);
  for (int i = 3; i < k; i++) {
    a += b;
    Decrement(bucket_start + (a & 7), (a >> 3) & 0x3f);
  }
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
void SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily, k>::
    Decrement(size_t group, int bit) {
#ifdef VERIFY_COUNT
    realCount[(group << 6) + bit]--;
#endif
    uint64_t m = data[group];
    uint64_t c = counts[group];
    if ((c & 0x8000000000000000ULL) != 0) {
        // an overflow entry
        size_t index = (size_t) (c & 0x0fffffffULL);
        size_t count = (size_t) (c >> 32) & 0x0fffffffULL;
        c -= 1ULL << 32;
        counts[group] = c;
        uint64_t n = overflow[index + bit / 8];
        overflow[index + bit / 8] = n - (1ULL << (bit * 8));
        n >>= 8 * (bit & 0xf);
        if ((n & 0xff) == 1) {
            data[group] &= ~(1ULL << bit);
        }
        if (count < 64) {
            // convert back to an inline entry, and free up the overflow entry
            uint64_t c2 = 0;
            for (int j = 63; j >= 0; j--) {
                int cj = (int) ((overflow[index + j / 8] >> (8 * j)) & 0xff);
                if (cj > 0) {
                    c2 = ((c2 << 1) | 1) << (cj - 1);
                }
            }
            counts[group] = c2;
            // free overflow
            overflow[index] = nextFreeOverflow;
            nextFreeOverflow = index;
        }
    } else {
        int bitsBefore = bitCount64(m & (0xffffffffffffffffULL >> (63 - bit)));
        int before = select64((c << 1) | 1, bitsBefore) - 1;
        int removeAt = max(0, before - 1);
        // remove the bit from the counter
        uint64_t mask = (1ULL << removeAt) - 1;
        uint64_t left = (c >> 1) & ~mask;
        uint64_t right= c & mask;
        counts[group] = left | right;
        uint64_t removed = (c >> removeAt) & 1;
        // possibly reset the data bit
        data[group] = m & ~(removed << bit);
    }
#ifdef VERIFY_COUNT
    VerifyCount(group, bit, __LINE__);
#endif
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
bool SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily, k>::
    Contain(const uint64_t key) const noexcept {
  const auto hash = hasher(key);
  const uint32_t bucket_start = reduce(rotl64(hash, 32), bucketCount) * 8;
  uint32_t a = (uint32_t)hash;
  char ok = 1;
  if (k >= 3) {
    ok &= data[bucket_start + ((a >> 0) & 7)] >> ((a >> 3) & 0x3f);
    ok &= data[bucket_start + ((a >> 9) & 7)] >> ((a >> 12) & 0x3f);
    ok &= data[bucket_start + ((a >> 18) & 7)] >> ((a >> 21) & 0x3f);
  }
  if (!ok) {
    return ok;
  }
  uint32_t b = (uint32_t)(hash >> 32);
  for (int i = 3; i < k; i++) {
    a += b;
    ok &= data[bucket_start + (a & 7)] >> ((a >> 3) & 63);
    if (!ok) {
      return ok;
    }
  }
  return ok;
}


// SuccinctCountingBlockedBloomRankFilter --------------------------------------------------------------------------------------

// #define VERIFY_COUNT

template <typename ItemType, size_t bits_per_item, typename HashFamily,
          int k = (int)((double)bits_per_item * 0.693147180559945 + 0.5)>
class SuccinctCountingBlockedBloomRankFilter {
private:
  const int bucketCount;
  HashFamily hasher;
  uint64_t *data;
  uint64_t *counts;
  uint64_t *overflow;
  size_t overflowLength;
  size_t nextFreeOverflow;
#ifdef VERIFY_COUNT
  uint8_t *realCount;
#endif

  void Increment(size_t group, int bit);
  void Decrement(size_t group, int bit);
  int ReadCount(size_t group, int bit);
#ifdef VERIFY_COUNT
  void VerifyCount(size_t group, int bit, int line);
#endif

public:
  explicit SuccinctCountingBlockedBloomRankFilter(const int capacity);
  ~SuccinctCountingBlockedBloomRankFilter() noexcept;
  void Add(const uint64_t key) noexcept;
  void Remove(const uint64_t key) noexcept;
  bool Contain(const uint64_t key) const noexcept;
  uint64_t SizeInBytes() const {
      return 2 * 64 * bucketCount + 8 * overflowLength;
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily, k>::
    SuccinctCountingBlockedBloomRankFilter(const int capacity)
    : bucketCount(capacity * bits_per_item / 512), hasher() {
  const size_t alloc_size = bucketCount * (512 / 8);
  const int malloc_failed =
      posix_memalign(reinterpret_cast<void **>(&data), 64, alloc_size);
  if (malloc_failed)
    throw ::std::bad_alloc();
  memset(data, 0, alloc_size);
  size_t arrayLength = bucketCount * 8;
  overflowLength = 100 + arrayLength / 100 * 36;
  counts = new uint64_t[arrayLength]();
  overflow = new uint64_t[overflowLength]();
#ifdef VERIFY_COUNT
  realCount = new uint8_t[arrayLength * 64]();
#endif
  nextFreeOverflow = 0;
  for (size_t i = 0; i < overflowLength; i += 8) {
      overflow[i] = i + 8;
  }
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily, k>::
    ~SuccinctCountingBlockedBloomRankFilter() noexcept {
  free(data);
  delete[] counts;
  delete[] overflow;
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
void SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily, k>::
    Add(const uint64_t key) noexcept {
  const auto hash = hasher(key);
  const uint32_t bucket_start = reduce(rotl64(hash, 32), bucketCount) * 8;
  uint32_t a = (uint32_t)hash;
  if (k >= 3) {
    Increment(bucket_start + ((a >> 0) & 7), (a >> 3) & 0x3f);
    Increment(bucket_start + ((a >> 9) & 7), (a >> 12) & 0x3f);
    Increment(bucket_start + ((a >> 18) & 7), (a >> 21) & 0x3f);
  }
  uint32_t b = (uint32_t)(hash >> 32);
  for (int i = 3; i < k; i++) {
    a += b;
    Increment(bucket_start + (a & 7), (a >> 3) & 0x3f);
  }
}

inline uint64_t getBit(int index) {
    return 1L << (index * 8);
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
void SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily, k>::
    Increment(size_t group, int x) {
#ifdef VERIFY_COUNT
    realCount[(group << 6) + x]++;
#endif
    uint64_t m = data[group];
    uint64_t c = counts[group];
    if ((c & 0x8000000000000000ULL) != 0) {
        // already an overflow
        size_t index = (int) (c & 0x0fffffff);
        c += 1L << 32;
        counts[group] = c;
        size_t bitIndex = x & 63;
        overflow[index + bitIndex / 8] += getBit(bitIndex);
        data[group] |= (1L << x);
#ifdef VERIFY_COUNT
        VerifyCount(group, x, __LINE__);
#endif
        return;
    }
    uint64_t d = (m >> x) & 1;
    if (d == 0 && c == 0) {
        data[group] |= 1L << x;
        return;
    }
    int bitsSet = bitCount64(m);
    int bitsBefore = x == 0 ? 0 : bitCount64(m << (64 - x));
    int insertAt = bitsBefore;
    if (d == 1) {
        int startLevel = 0;
        uint64_t bitsForLevel;
        while (true) {
            uint64_t levelMask = ((1L << bitsSet) - 1) << startLevel;
            bitsForLevel = c & levelMask;
            if (((c >> insertAt) & 1) == 0) {
                break;
            }
            startLevel += bitsSet;
            if (startLevel >= 64) {
                // convert to overflow later
                insertAt = 64;
                break;
            }
            bitsSet = bitCount64(bitsForLevel);
            bitsBefore = insertAt == 0 ? 0 : bitCount64(bitsForLevel << (64 - insertAt));
            insertAt = startLevel + bitsBefore;
        }
        // bit is not set: set it, and insert a space in the next level if needed
        c |= 1L << insertAt;
        int bitsBeforeLevel = insertAt == 0 ? 0 : bitCount64(bitsForLevel << (64 - insertAt));
        int bitsSetLevel = bitCount64(bitsForLevel);
        insertAt = startLevel + bitsSet + bitsBeforeLevel;
        bitsSet = bitsSetLevel;
    }
    // insert a space
    uint64_t mask = (1L << insertAt) - 1;
    uint64_t left = c & ~mask;
    uint64_t right = c & mask;
    c = (left << 1) | right;
    if (insertAt >= 64 || (c & 0x8000000000000000L) != 0) {
        // an overflow entry, or overflowing now
        int index = nextFreeOverflow;
        nextFreeOverflow = (int) overflow[index];
        for (int i = 0; i < 8; i++) {
            overflow[index + i] = 0;
        }
        // convert to a pointer
        uint64_t count = 1;
        for (int i = 0; i < 64; i++) {
            int n = ReadCount(group, i);
            count += n;
            overflow[index + i / 8] += n * getBit(i);
        }
        c = 0x8000000000000000L | (count << 32) | index;
        data[group] |= 1L << x;
        counts[group] = c;
        int bitIndex = x & 63;
        overflow[index + bitIndex / 8] += getBit(bitIndex);
#ifdef VERIFY_COUNT
        VerifyCount(group, x, __LINE__);
#endif
        return;
    }
    data[group] |= 1L << x;
    counts[group] = c;
#ifdef VERIFY_COUNT
    VerifyCount(group, x, __LINE__);
#endif
}

#ifdef VERIFY_COUNT
template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
void SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily, k>::
    VerifyCount(size_t group, int bit, int line) {
    for(int b = 0; b < 64; b++) {
        if (realCount[(group << 6) + b] != ReadCount(group, b)) {
            ::std::cout << "group " << group << " bit " << b <<
            " expected " << (int) realCount[(group << 6) + b] <<
            " got " << ReadCount(group, b) <<
            " at bit " << bit << " at line " << line << "\n";
        }
    }
}
#endif

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
int SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily, k>::
    ReadCount(size_t group, int x) {
    uint64_t m = data[group];
    uint64_t d = (m >> x) & 1;
    if (d == 0) {
        return 0;
    }
    uint64_t c = counts[group];
    if ((c & 0x8000000000000000L) != 0) {
        int index = (int) (c & 0x0fffffff);
        int bitIndex = x & 63;
        uint64_t n = overflow[index + bitIndex / 8];
        n >>= 8 * (bitIndex & 7);
        return (int) (n & 0xff);
    }
    if (c == 0) {
        return 1;
    }
    int bitsSet = bitCount64(m);
    x &= 63;
    int bitsBefore = x == 0 ? 0 : bitCount64(m << (64 - x));
    int count = 1;
    int insertAt = bitsBefore;
    while (true) {
        // the mask for the current level
        uint64_t levelMask = (1L << bitsSet) - 1;
        // the relevant bits for the current level
        uint64_t bitsForLevel = c & levelMask;
        if (((c >> insertAt) & 1) == 1) {
            count++;
            // at this level, the bit is already set: loop until it's not set
            c >>= bitsSet;
            bitsSet = bitCount64(bitsForLevel);
            bitsBefore = insertAt == 0 ? 0 : bitCount64(bitsForLevel << (64 - insertAt));
            insertAt = bitsBefore;
        } else {
            break;
        }
        if (count > 16) {
            // unexpected
            ::std::cout << "group- " << group << " count " << count << "\n";
        }
    }
    return count;
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
void SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily, k>::
    Remove(const uint64_t key) noexcept {
  const auto hash = hasher(key);
  const uint32_t bucket_start = reduce(rotl64(hash, 32), bucketCount) * 8;
  uint32_t a = (uint32_t)hash;
  if (k >= 3) {
    Decrement(bucket_start + ((a >> 0) & 7), (a >> 3) & 0x3f);
    Decrement(bucket_start + ((a >> 9) & 7), (a >> 12) & 0x3f);
    Decrement(bucket_start + ((a >> 18) & 7), (a >> 21) & 0x3f);
  }
  uint32_t b = (uint32_t)(hash >> 32);
  for (int i = 3; i < k; i++) {
    a += b;
    Decrement(bucket_start + (a & 7), (a >> 3) & 0x3f);
  }
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
void SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily, k>::
    Decrement(size_t group, int x) {
#ifdef VERIFY_COUNT
    realCount[(group << 6) + x]--;
#endif
    uint64_t m = data[group];
    uint64_t c = counts[group];
    if ((c & 0x8000000000000000L) != 0) {
        // an overflow entry
        int count = (int) (c >> 32) & 0x0fffffff;
        c -= 1L << 32;
        counts[group] = c;
        int index = (int) (c & 0x0fffffff);
        int bitIndex = x & 63;
        uint64_t n = overflow[index + bitIndex / 8];
        overflow[index + bitIndex / 8] = n - getBit(bitIndex);
        n >>= 8 * (bitIndex & 7);
        if ((n & 0xff) == 1) {
            data[group] &= ~(1L << x);
        }
        if (count < 64) {
            // convert back to an inline entry, and free up the overflow entry
            int temp[64];
            int count2 = 0;
            for(int j = 0; j < 64; j++) {
                int cj = (int) ((overflow[index + j / 8] >> (8 * (j & 7))) & 0xff);
                count2 += cj;
                temp[j] = cj;
            }
            uint64_t c2 = 0;
            int off = 0;
            while (count2 > 0) {
                for (int i = 0; i < 64; i++) {
                    int t = temp[i];
                    if (t > 0) {
                        temp[i]--;
                        count2--;
                        c2 |= ((t > 1) ? 1L : 0L) << off;
                        off++;
                    }
                }
            }
            counts[group] = c2;
            // freeOverflow(index);
            overflow[index] = nextFreeOverflow;
            nextFreeOverflow = index;
#ifdef VERIFY_COUNT
            VerifyCount(group, x, __LINE__);
#endif
        }
        return;
    }
    // number of bits in the counter at this level
    int bitsSet = bitCount64(m);
    // number of bits before the bit to test (at the current level)
    int bitsBefore = x == 0 ? 0 : bitCount64(m << (64 - x));
    // the point where the bit should be removed (remove 0), or reset
    int removeAt = bitsBefore;
    uint64_t d = (c >> bitsBefore) & 1;
    if (d == 1) {
        // bit is set: loop until it's not set
        int startLevel = 0;
        uint64_t bitsForLevel;
        int resetAt = removeAt;
        while (true) {
            // the mask for the current level
            uint64_t levelMask = ((1L << bitsSet) - 1) << startLevel;
            // the relevant bits for the current level
            bitsForLevel = c & levelMask;
            if (((c >> removeAt) & 1) == 0) {
                break;
            }
            // at this level, the bit is already set: loop until it's not set
            startLevel += bitsSet;
            bitsSet = bitCount64(bitsForLevel);
            bitsBefore = removeAt == 0 ? 0 : bitCount64(bitsForLevel << (64 - removeAt));
            resetAt = removeAt;
            removeAt = startLevel + bitsBefore;
            if (removeAt > 63) {
                break;
            }
        }
        c ^= 1L << resetAt;
    }
    if (removeAt < 64) {
        // remove the bit from the counter
        uint64_t mask = (1L << removeAt) - 1;
        uint64_t left = (c >> 1) & ~mask;
        uint64_t right= c & mask;
        c = left | right;
    }
    counts[group] = c;
    // possibly reset the data bit
    data[group] = m & ~((d==0?1L:0L) << x);
#ifdef VERIFY_COUNT
    VerifyCount(group, x, __LINE__);
#endif
}

template <typename ItemType, size_t bits_per_item, typename HashFamily, int k>
bool SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily, k>::
    Contain(const uint64_t key) const noexcept {
  const auto hash = hasher(key);
  const uint32_t bucket_start = reduce(rotl64(hash, 32), bucketCount) * 8;
  uint32_t a = (uint32_t)hash;
  char ok = 1;
  if (k >= 3) {
    ok &= data[bucket_start + ((a >> 0) & 7)] >> ((a >> 3) & 0x3f);
    ok &= data[bucket_start + ((a >> 9) & 7)] >> ((a >> 12) & 0x3f);
    ok &= data[bucket_start + ((a >> 18) & 7)] >> ((a >> 21) & 0x3f);
  }
  if (!ok) {
    return ok;
  }
  uint32_t b = (uint32_t)(hash >> 32);
  for (int i = 3; i < k; i++) {
    a += b;
    ok &= data[bucket_start + (a & 7)] >> ((a >> 3) & 63);
    if (!ok) {
      return ok;
    }
  }
  return ok;
}

}
#endif