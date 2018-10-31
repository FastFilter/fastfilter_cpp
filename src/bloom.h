#ifndef BLOOM_FILTER_BLOOM_FILTER_H_
#define BLOOM_FILTER_BLOOM_FILTER_H_

#include <assert.h>
#include <algorithm>

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
    return (uint32_t) (((uint64_t) hash * n) >> 32);
}

static size_t getBestK(size_t bitsPerItem) {
    return max(1, (int) round((double) bitsPerItem * log(2)));
}

inline uint64_t getBit(uint32_t index) {
    return 1L << (index & 63);
}

template <typename ItemType, size_t bits_per_item,
    typename HashFamily = TwoIndependentMultiplyShift,
    int k = (int) ((double) bits_per_item * 0.693147180559945 + 0.5)>
class BloomFilter {

  uint64_t *data;
  size_t size;
  size_t arrayLength;
  size_t bitCount;
  int kk;
  HashFamily hasher;

  double BitsPerItem() const { return k; }

 public:
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
  Status AddAll(const vector<ItemType> data, const size_t start, const size_t end);

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

template <typename ItemType, size_t bits_per_item,
    typename HashFamily, int k>
Status BloomFilter<ItemType, bits_per_item, HashFamily, k>::Add(
    const ItemType &key) {
    uint64_t hash = hasher(key);
    uint32_t a = (uint32_t) (hash >> 32);
    uint32_t b = (uint32_t) hash;
    for (int i = 0; i < k; i++) {
        // int index = reduce(a, this->bitCount);
        // data[index >> 6] |= getBit(index);
        // reworked to avoid overflows
        // use the fact that reduce is not very sensitive to lower bits of a
        data[reduce(a, this->arrayLength)] |= getBit(a);
        a += b;
    }
    return Ok;
}

const int blockShift = 15;
const int blockLen = 1 << blockShift;

void applyBlock(uint32_t* tmp, int block, int len, uint64_t *data) {
    for (int i = 0; i < len; i++) {
        uint32_t index = tmp[(block << blockShift) + i];
        data[index >> 6] |= getBit(index);
    }
}

template <typename ItemType, size_t bits_per_item,
    typename HashFamily, int k>
Status BloomFilter<ItemType, bits_per_item, HashFamily, k>::AddAll(
    const vector<ItemType> keys, const size_t start, const size_t end) {
    int blocks = 1 + arrayLength / blockLen;
    uint32_t* tmp = new uint32_t[blocks * blockLen];
    int* tmpLen = new int[blocks]();
    for(size_t i = start; i < end; i++) {
        uint64_t key = keys[i];
        uint64_t hash = hasher(key);
        uint32_t a = (uint32_t) (hash >> 32);
        uint32_t b = (uint32_t) hash;
        for (int j = 0; j < k; j++) {
            int index = reduce(a, this->arrayLength);
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

template <typename ItemType, size_t bits_per_item,
    typename HashFamily, int k>
Status BloomFilter<ItemType, bits_per_item, HashFamily, k>::Contain(
    const ItemType &key) const {
    uint64_t hash = hasher(key);
    uint32_t a = (uint32_t) (hash >> 32);
    uint32_t b = (uint32_t) hash;
    for (int i = 0; i < k; i++) {
        // int index = reduce(a, this->bitCount);
        // if ((data[index >> 6] & getBit(index)) == 0) {
        //     return NotFound;
        // }
        // reworked to avoid overflows
        if ((data[reduce(a, this->arrayLength)] & getBit(a)) == 0) {
            return NotFound;
        }
        a += b;
    }
    return Ok;
}

template <typename ItemType, size_t bits_per_item,
    typename HashFamily, int k>
std::string BloomFilter<ItemType, bits_per_item, HashFamily, k>::Info() const {
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

template<size_t blocksize, int k, typename HashFamily = ::hashing::TwoIndependentMultiplyShift>
class SimpleBlockFilter {
 private:
  // The filter is divided up into Buckets:
  using Bucket = uint64_t[blocksize];

  const int bucketCount;

  Bucket* directory_;

  HashFamily hasher_;

 public:
  // Consumes at most (1 << log_heap_space) bytes on the heap:
  explicit SimpleBlockFilter(const int bits);
  ~SimpleBlockFilter() noexcept;
  void Add(const uint64_t key) noexcept;

  bool Find(const uint64_t key) const noexcept;
  uint64_t SizeInBytes() const { return sizeof(Bucket) * bucketCount; }


};

template<size_t blocksize,int k,  typename HashFamily>
SimpleBlockFilter<blocksize,k,HashFamily>::SimpleBlockFilter(const int capacity)
  : bucketCount(capacity * k / (  blocksize * 64)),
    directory_(nullptr),
    hasher_() {
  const size_t alloc_size = bucketCount * sizeof(Bucket);
  const int malloc_failed =
      posix_memalign(reinterpret_cast<void**>(&directory_), 64, alloc_size);
  if (malloc_failed) throw ::std::bad_alloc();
  memset(directory_, 0, alloc_size);
}

template<size_t blocksize, int k, typename HashFamily>
SimpleBlockFilter<blocksize,k,HashFamily>::~SimpleBlockFilter() noexcept {
  free(directory_);
  directory_ = nullptr;
}

static inline uint64_t rotl64(uint64_t n, unsigned int c) {
    // assumes width is a power of 2
    const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
    // assert ( (c<=mask) &&"rotate by type width or more");
    c &= mask;
    return (n << c) | ( n >> ((-c) & mask));
}

template<size_t blocksize,int k,  typename HashFamily>
inline void
SimpleBlockFilter<blocksize,k, HashFamily>::Add(const uint64_t key) noexcept {
  const auto hash = hasher_(key);
  const uint32_t bucket_idx = reduce(rotl64(hash, 32), bucketCount);
  Bucket * bucket = directory_ + bucket_idx;
  uint32_t a = (uint32_t) (hash >> 32);
  uint32_t b = (uint32_t) hash;
  for (int i = 0; i < k; i++) {
        ((uint64_t *)bucket)[reduce(a, blocksize)] |= getBit(a);
        a += b;
  }

}


template<size_t blocksize, int k, typename HashFamily>
inline bool
SimpleBlockFilter<blocksize,k,HashFamily>::Find(const uint64_t key) const noexcept {
  const auto hash = hasher_(key);
  const uint32_t bucket_idx = reduce(rotl64(hash, 32), bucketCount);
  const Bucket * bucket = directory_ + bucket_idx;
  uint32_t a = (uint32_t) (hash >> 32);
  uint32_t b = (uint32_t) hash;
  for (int i = 0; i < k; i++) {
        if ((((const uint64_t *)bucket)[reduce(a, blocksize)] & getBit(a)) == 0) {
            return false;
        }
        a += b;
  }
  return true;
}



}  // namespace bloomfilter
#endif  // BLOOM_FILTER_BLOOM_FILTER_H_
