#ifndef BLOOM_FILTER_BLOOM_FILTER_H_
#define BLOOM_FILTER_BLOOM_FILTER_H_

#include <assert.h>
#include <algorithm>

#include "debug.h"
#include "hashutil.h"
#include "printutil.h"

using namespace std;
using namespace cuckoofilter;

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

// https://stackoverflow.com/questions/30052316/find-next-prime-number-algorithm
bool isPrime(int number) {
    if (number == 2 || number == 3)
        return true;
    if (number % 2 == 0 || number % 3 == 0)
        return false;
    int divisor = 6;
    while (divisor * divisor - 2 * divisor + 1 <= number) {
        if (number % (divisor - 1) == 0)
            return false;
        if (number % (divisor + 1) == 0)
            return false;
        divisor += 6;
    }
    return true;
}

int nextPrime(int a) {
    while (!isPrime(++a));
    return a;
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
    this->bitCount = nextPrime(n * bits_per_item);
    this->arrayLength = (bitCount + 63) / 64;
    data = new uint64_t[arrayLength];
    std::fill_n(data, arrayLength, 0);
  }

  ~BloomFilter() { delete[] data; }

  // Add an item to the filter.
  Status Add(const ItemType &item);

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
}  // namespace bloomfilter
#endif  // BLOOM_FILTER_BLOOM_FILTER_H_
