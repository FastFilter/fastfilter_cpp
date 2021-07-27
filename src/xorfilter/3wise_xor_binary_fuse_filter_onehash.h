#ifndef THREEWISE_XOR_BINARY_FUSE_FILTER_XOR_FILTER_ONEHASH_H_
#define THREEWISE_XOR_BINARY_FUSE_FILTER_XOR_FILTER_ONEHASH_H_
#include "xor_binary_fuse_filter.h"
namespace xorbinaryfusefilter_onehash {
// status returned by a xor filter operation
enum Status {
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

inline uint64_t rotl64(uint64_t n, unsigned int c) {
  // assumes width is a power of 2
  const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
  // assert ( (c<=mask) &&"rotate by type width or more");
  c &= mask;
  return (n << c) | (n >> ((-c) & mask));
}

__attribute__((always_inline)) inline uint32_t reduce(uint32_t hash,
                                                      uint32_t n) {
  // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
  return (uint32_t)(((uint64_t)hash * n) >> 32);
}

template <typename ItemType, typename FingerprintType,
          typename HashFamily = SimpleMixSplit>
class XorBinaryFuseFilter {
public:
  size_t size;
  size_t arrayLength;
  size_t segmentCount;
  size_t segmentCountLength;
  size_t segmentLength;
  size_t segmentLengthMask;
  static constexpr size_t arity = 3;
  FingerprintType *fingerprints;

  HashFamily *hasher;

  inline FingerprintType fingerprint(const uint64_t hash) const {
    return (FingerprintType)hash;
  }

  inline __attribute__((always_inline)) size_t getHashFromHash(uint64_t hash,
                                                               int index) {
    __uint128_t x = (__uint128_t)hash * (__uint128_t)segmentCountLength;
    uint64_t h = (uint64_t)(x >> 64);
    h += index * segmentLength;
    // keep the lower 36 bits
    uint64_t hh = hash & ((1UL << 36) - 1);
    // index 0: right shift by 36; index 1: right shift by 18; index 2: no shift
    h ^= (size_t)((hh >> (36 - 18 * index)) & segmentLengthMask);
    return h;
  }

  explicit XorBinaryFuseFilter(const size_t size) {
    hasher = new HashFamily();
    this->size = size;
    this->segmentLength = calculateSegmentLength(arity, size);
    // the current implementation hardcodes a 18-bit limit to
    // to the segment length.
    if (this->segmentLength > (1 << 18)) {
      this->segmentLength = (1 << 18);
    }
    double sizeFactor = calculateSizeFactor(arity, size);
    size_t capacity = size * sizeFactor;
    size_t segmentCount =
        (capacity + segmentLength - 1) / segmentLength - (arity - 1);
    this->arrayLength = (segmentCount + arity - 1) * segmentLength;
    this->segmentLengthMask = this->segmentLength - 1;
    this->segmentCount =
        (this->arrayLength + this->segmentLength - 1) / this->segmentLength;
    this->segmentCount =
        this->segmentCount <= arity - 1 ? 1 : this->segmentCount - (arity - 1);
    this->arrayLength = (this->segmentCount + arity - 1) * this->segmentLength;
    this->segmentCountLength = this->segmentCount * this->segmentLength;
    fingerprints = new FingerprintType[arrayLength]();
    std::fill_n(fingerprints, arrayLength, 0);
  }

  ~XorBinaryFuseFilter() {
    delete[] fingerprints;
    delete hasher;
  }

  Status AddAll(const vector<ItemType> &data, const size_t start,
                const size_t end) {
    return AddAll(data.data(), start, end);
  }

  Status AddAll(const ItemType *data, const size_t start, const size_t end);

  // Report if the item is inserted, with false positive rate.
  Status Contain(const ItemType &item) const;

  /* methods for providing stats  */
  // summary infomation
  std::string Info() const;

  // number of current inserted items;
  size_t Size() const { return size; }

  // size of the filter in bytes.
  size_t SizeInBytes() const { return arrayLength * sizeof(FingerprintType); }
};

struct t2val {
  uint8_t count;
  uint32_t indexA;
  uint32_t indexB;
  uint16_t fingerprint;
};

typedef struct t2val t2val_t;

struct reverse {
  uint32_t index;
  uint32_t indexA;
  uint32_t indexB;
  uint16_t fingerprint;
};

typedef struct reverse reverse_t;


template <typename ItemType, typename FingerprintType, typename HashFamily>
Status XorBinaryFuseFilter<ItemType, FingerprintType, HashFamily>::AddAll(
    const ItemType *keys, const size_t start, const size_t end) {

  uint64_t *partialSort = new uint64_t[size + 1];
  partialSort[size] = 1;
  
  reverse_t *reverseOrder = new reverse_t[size];
  
  size_t reverseOrderPos;

  t2val_t *t2vals = new t2val_t[arrayLength];

  size_t *alone = new size_t[arrayLength];
  size_t hashIndex{0};

  while (true) {
    memset(t2vals, 0, sizeof(t2val_t[arrayLength]));

    // counting sort

    memset(partialSort, 0, sizeof(uint64_t[size]));

    int blockBits = 1;
    while((size_t(1)<<blockBits) < segmentCount) { blockBits++; }
    size_t block = size_t(1) << blockBits;
    size_t *startPos = new size_t[block];
    for(uint32_t i = 0; i < uint32_t(1) << blockBits; i++) { startPos[i] = i * size / block; }
    
    for (size_t i = start; i < end; i++) {
      uint64_t k = keys[i];
      uint64_t hash = (*hasher)(k);
      size_t segment_index = hash >> (64 - blockBits);
      // We only overwrite when the hash was zero. Zero hash values
      // may be misplaced (unlikely).
      while(partialSort[startPos[segment_index]] != 0) {
        segment_index++;
        segment_index &= (size_t(1) << blockBits) - 1;
      }
      partialSort[startPos[segment_index]] = hash;
      startPos[segment_index]++;
    }
    
    for (size_t i = 0; i < size; i++) {
      uint64_t hash = partialSort[i];
      int index0 = getHashFromHash(hash, 0);
      int index1 = getHashFromHash(hash, 1);
      int index2 = getHashFromHash(hash, 2);
      FingerprintType fp = fingerprint(hash);
      t2vals[index0].indexA ^= index1;
      t2vals[index0].indexB ^= index2;
      t2vals[index1].indexA ^= index2;
      t2vals[index1].indexB ^= index0;
      t2vals[index2].indexA ^= index0;
      t2vals[index2].indexB ^= index1;
      t2vals[index0].fingerprint ^= fp;
      t2vals[index1].fingerprint ^= fp;
      t2vals[index2].fingerprint ^= fp;
      t2vals[index0].count++;
      t2vals[index1].count++;
      t2vals[index2].count++;
    }
    delete[] startPos;

    reverseOrderPos = 0;
    size_t alonePos = 0;
    for (size_t i = 0; i < arrayLength; i++) {
      int if1 = t2vals[i].count == 1 ? 1 : 0;
      alone[alonePos] = i;
      alonePos += if1;
    }

    while (alonePos > 0) {
      alonePos--;
      size_t index = alone[alonePos];
      if (t2vals[index].count == 1) {
        size_t indexA = t2vals[index].indexA;
        size_t indexB = t2vals[index].indexB;
        uint32_t fp = t2vals[index].fingerprint;
        reverseOrder[reverseOrderPos].index = index;
        reverseOrder[reverseOrderPos].indexA = indexA;
        reverseOrder[reverseOrderPos].indexB = indexB;
        reverseOrder[reverseOrderPos].fingerprint = fp;
        reverseOrderPos++;
        t2vals[indexA].count--;
        t2vals[indexA].indexA ^= indexB;
        t2vals[indexA].indexB ^= index;
        t2vals[indexA].fingerprint ^= fp;
        t2vals[indexB].count--;
        t2vals[indexB].indexA ^= index;
        t2vals[indexB].indexB ^= indexA;
        t2vals[indexB].fingerprint ^= fp;
        int ifA = t2vals[indexA].count == 1 ? 1 : 0;
        alone[alonePos] = indexA;
        alonePos += ifA;
        int ifB = t2vals[indexB].count == 1 ? 1 : 0;
        alone[alonePos] = indexB;
        alonePos += ifB;
      }
    }

    if (reverseOrderPos == size) {
      break;
    }
    hashIndex++;
    // use a new random numbers
    delete hasher;
    hasher = new HashFamily();
  }
  delete[] alone;
  delete[] t2vals;

  for (int i = reverseOrderPos - 1; i >= 0; i--) {
    reverse_t next = reverseOrder[i];
    FingerprintType fp = next.fingerprint;
    size_t index = next.index;
    size_t indexA = next.indexA;
    size_t indexB = next.indexB;
    fingerprints[index] = fingerprints[indexA] ^ fingerprints[indexB] ^ fp;
  }
  delete[] reverseOrder;

  return Ok;
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
Status XorBinaryFuseFilter<ItemType, FingerprintType, HashFamily>::Contain(
    const ItemType &key) const {
  uint64_t hash = (*hasher)(key);
  FingerprintType f = fingerprint(hash);
  __uint128_t x = (__uint128_t)hash * (__uint128_t)segmentCountLength;
  int h0 = (uint64_t)(x >> 64);
  int h1 = h0 + segmentLength;
  int h2 = h1 + segmentLength;
  uint64_t hh = hash;
  h1 ^= (size_t)((hh >> 18) & segmentLengthMask);
  h2 ^= (size_t)((hh)&segmentLengthMask);
  f ^= fingerprints[h0] ^ fingerprints[h1] ^ fingerprints[h2];
  return f == 0 ? Ok : NotFound;
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
std::string
XorBinaryFuseFilter<ItemType, FingerprintType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "XorBinaryFuseFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  return ss.str();
}
} // namespace xorbinaryfusefilter_onehash
#endif