#ifndef BLOOM_FILTER_BLOOM_FILTER_H_
#define BLOOM_FILTER_BLOOM_FILTER_H_

#include <algorithm>
#include <assert.h>
#include <sstream>

#include "hashutil.h"

using namespace std;
using namespace hashing;

namespace bloomfilter {
// status returned by a Bloom filter operation
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

/**
* Given a value "word", produces an integer in [0,p) without division.
* The function is as fair as possible in the sense that if you iterate
* through all possible values of "word", then you will generate all
* possible outputs as uniformly as possible.
*/
static inline uint32_t fastrange32(uint32_t word, uint32_t p) {
  // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
	return (uint32_t)(((uint64_t)word * (uint64_t)p) >> 32);
}

#if defined(_MSC_VER) && defined (_WIN64)
#include <intrin.h>// should be part of all recent Visual Studio
#pragma intrinsic(_umul128)
#endif // defined(_MSC_VER) && defined (_WIN64)


/**
* Given a value "word", produces an integer in [0,p) without division.
* The function is as fair as possible in the sense that if you iterate
* through all possible values of "word", then you will generate all
* possible outputs as uniformly as possible.
*/
static inline uint64_t fastrange64(uint64_t word, uint64_t p) {
  // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
#ifdef __SIZEOF_INT128__ // then we know we have a 128-bit int
	return (uint64_t)(((__uint128_t)word * (__uint128_t)p) >> 64);
#elif defined(_MSC_VER) && defined(_WIN64)
	// supported in Visual Studio 2005 and better
	uint64_t highProduct;
	_umul128(word, p, &highProduct); // ignore output
	return highProduct;
	unsigned __int64 _umul128(
		unsigned __int64 Multiplier,
		unsigned __int64 Multiplicand,
		unsigned __int64 *HighProduct
	);
#else
	return word % p; // fallback
#endif // __SIZEOF_INT128__
}


#ifndef UINT32_MAX
#define UINT32_MAX  (0xffffffff)
#endif // UINT32_MAX

/**
* Given a value "word", produces an integer in [0,p) without division.
* The function is as fair as possible in the sense that if you iterate
* through all possible values of "word", then you will generate all
* possible outputs as uniformly as possible.
*/
static inline size_t fastrangesize(uint64_t word, size_t p) {
#if (SIZE_MAX == UINT32_MAX)
	return (size_t)fastrange32(word, p);
#else // assume 64-bit
	return (size_t)fastrange64(word, p);
#endif // SIZE_MAX == UINT32_MAX
}

static inline size_t getBestK(size_t bitsPerItem) {
  return max(1, (int)round((double)bitsPerItem * log(2)));
}

inline uint64_t getBit(uint32_t index) { return 1L << (index & 63); }

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily = SimpleMixSplit,
          int k = (int)((double)bits_per_item * 0.693147180559945 + 0.5)>
class BloomFilter {
public:

  uint64_t *data;
  size_t size;
  size_t arrayLength;
  size_t bitCount;
  int kk;
  HashFamily hasher;

  double BitsPerItem() const { return k; }

  explicit BloomFilter(const size_t n) : hasher() {
    this->size = 0;
    this->kk = getBestK(bits_per_item);
    this->bitCount = n * bits_per_item;
    this->arrayLength = (bitCount + 63) / 64;
    data = new uint64_t[arrayLength];
    std::fill_n(data, arrayLength, 0);
  }

  ~BloomFilter() { delete[] data; }

  // Add an item to the filter.
  Status Add(const ItemType &item);

  // Add multiple items to the filter.
  Status AddAll(const vector<ItemType>& data, const size_t start,
                const size_t end) {
    return AddAll(data.data(),start,end);

  }
  Status AddAll(const ItemType* data, const size_t start,
                const size_t end);
  // Report if the item is inserted, with false positive rate.
  Status Contain(const ItemType &item) const;

  /* methods for providing stats  */
  // summary infomation
  std::string Info() const;

  // number of current inserted items;
  size_t Size() const { return size; }

  // size of the filter in bytes.
  size_t SizeInBytes() const { return arrayLength * 8; }
};

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
Status BloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::Add(
    const ItemType &key) {
  uint64_t hash = hasher(key);
  uint64_t a = (hash >> 32) | (hash << 32);
  uint64_t b = hash;
  for (int i = 0; i < k; i++) {
    // int index = reduce(a, this->bitCount);
    // data[index >> 6] |= getBit(index);
    // reworked to avoid overflows
    // use the fact that reduce is not very sensitive to lower bits of a
    data[fastrangesize(a, this->arrayLength)] |= getBit(a);
    a += b;
  }
  return Ok;
}

const int blockShift = 15;
const int blockLen = 1 << blockShift;

void applyBlock(uint32_t *tmp, int block, int len, uint64_t *data) {
  for (int i = 0; i < len; i++) {
    uint32_t index = tmp[(block << blockShift) + i];
    data[index >> 6] |= getBit(index);
  }
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
Status BloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::AddAll(
    const ItemType* keys, const size_t start, const size_t end) {
  // we have that AddAll assumes that arrayLength << 6 is a
  // 32-bit integer
  if(arrayLength > 0x3ffffff) {
    for(size_t i = start; i < end; i++) {
      Add(keys[i]);
    }
    return Ok;
  }
  int blocks = 1 + arrayLength / blockLen;
  uint32_t *tmp = new uint32_t[blocks * blockLen];
  int *tmpLen = new int[blocks]();
  for (size_t i = start; i < end; i++) {
    uint64_t key = keys[i];
    uint64_t hash = hasher(key);
    uint64_t a = (hash >> 32) | (hash << 32);
    uint64_t b = hash;
    for (int j = 0; j < k; j++) {
      int index = fastrangesize(a, this->arrayLength);
      int block = index >> blockShift;
      int len = tmpLen[block];
      tmp[(block << blockShift) + len] = (index << 6) + (a & 63);
      tmpLen[block] = len + 1;
      if (len + 1 == blockLen) {
        applyBlock(tmp, block, len + 1, data);
        tmpLen[block] = 0;
      }
      a += b;
    }
  }
  for (int block = 0; block < blocks; block++) {
    applyBlock(tmp, block, tmpLen[block], data);
  }
  delete[] tmp;
  delete[] tmpLen;
  return Ok;
}

char bittest64(const uint64_t *t, uint64_t bit) {
  return (*t & (1L << (bit & 63))) != 0;
}
template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
Status BloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::Contain(
    const ItemType &key) const {
  uint64_t hash = hasher(key);
  uint64_t a = (hash >> 32) | (hash << 32);
  uint64_t b = hash;
  if (branchless && k >= 3) {
    int b0 = data[fastrangesize(a, this->arrayLength)] >> (a & 63);
    a += b;
    int b1 = data[fastrangesize(a, this->arrayLength)] >> (a & 63);
    a += b;
    int b2 = data[fastrangesize(a, this->arrayLength)] >> (a & 63);
    if ((b0 & b1 & b2 & 1) == 0) {
        return NotFound;
    }
    for (int i = 3; i < k; i++) {
      a += b;
      if (((data[fastrangesize(a, this->arrayLength)] >> (a & 63)) & 1) == 0) {
          return NotFound;
      }
    }
    return Ok;
  }
  for (int i = 0; i < k; i++) {
    if ((data[fastrangesize(a, this->arrayLength)] & getBit(a)) == 0) {
      return NotFound;
    }
    a += b;
  }
  return Ok;
}

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily, int k>
std::string
BloomFilter<ItemType, bits_per_item, branchless, HashFamily, k>::Info() const {
  std::stringstream ss;
  ss << "BloomFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  if (Size() > 0) {
    ss << "\t\tk:   " << BitsPerItem() << "\n";
  } else {
    ss << "\t\tk:   N/A\n";
  }
  return ss.str();
}



/***************
 * Simple block filter (naive implementation)
 ***************/

template <size_t blocksize, int k,
          typename HashFamily = ::hashing::SimpleMixSplit>
class SimpleBlockFilter {
private:
  const size_t arrayLength;
  uint64_t* data;
  HashFamily hasher_;
public:
  // Consumes at most (1 << log_heap_space) bytes on the heap:
  explicit SimpleBlockFilter(const int bits);
  ~SimpleBlockFilter() noexcept;
  void Add(const uint64_t key) noexcept;
  bool Find(const uint64_t key) const noexcept;
  uint64_t SizeInBytes() const {
    return arrayLength * 8;
  }
};

template <size_t blocksize, int k, typename HashFamily>
SimpleBlockFilter<blocksize, k, HashFamily>::SimpleBlockFilter(
    const int capacity)
    : arrayLength((capacity * 10) / 64 + 8),
      hasher_() {
  data = new uint64_t[arrayLength]();
}

template <size_t blocksize, int k, typename HashFamily>
SimpleBlockFilter<blocksize, k, HashFamily>::~SimpleBlockFilter() noexcept {
  delete[] data;
  data = nullptr;
}

static inline uint64_t rotl64(uint64_t n, unsigned int c) {
  // assumes width is a power of 2
  const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
  c &= mask;
  return (n << c) | (n >> ((-c) & mask));
}

template <size_t blocksize, int k, typename HashFamily>
inline void
SimpleBlockFilter<blocksize, k, HashFamily>::Add(const uint64_t key) noexcept {
  const auto hash = hasher_(key);
  const uint32_t idx = reduce(hash, arrayLength);
  uint64_t *bucket = data + idx;
  uint64_t m1 = 1L << hash;
  uint64_t m2 = 1L << (hash >> 8);
  uint64_t m = m1 | m2;
  *bucket |= m;

}
template <size_t blocksize, int k, typename HashFamily>
inline bool
SimpleBlockFilter<blocksize, k, HashFamily>::Find(const uint64_t key) const
    noexcept {
  const auto hash = hasher_(key);
  const uint32_t idx = reduce(hash, arrayLength);
  uint64_t *bucket = data + idx;
  uint64_t m1 = 1L << hash;
  uint64_t m2 = 1L << (hash >> 8);
  uint64_t m = m1 | m2;
  return !((m & *bucket) - m);
}

} // namespace bloomfilter
#endif // BLOOM_FILTER_BLOOM_FILTER_H_
