#ifndef THREEWISE_XOR_BINARY_FUSE_FILTER_XOR_FILTER_LOWMEM_H_
#define THREEWISE_XOR_BINARY_FUSE_FILTER_XOR_FILTER_LOWMEM_H_
#include "xor_binary_fuse_filter.h"

/**
 * As of July 2021, the lowmem versions of the binary fuse filters are
 * the recommended defaults.
 */
namespace xorbinaryfusefilter_lowmem {
// status returned by a xor filter operation
enum Status {
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

__attribute__((always_inline)) inline uint32_t reduce(uint32_t hash,
                                                      uint32_t n) {
  // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
  return (uint32_t)(((uint64_t)hash * n) >> 32);
}

__attribute__((always_inline)) inline uint8_t mod3(uint8_t x) {
    if (x > 2) {
        x -= 3;
    }
    return x;
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
  size_t hashIndex{0};

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

template <typename ItemType, typename FingerprintType, typename HashFamily>
Status XorBinaryFuseFilter<ItemType, FingerprintType, HashFamily>::AddAll(
    const ItemType *keys, const size_t start, const size_t end) {

  uint64_t *reverseOrder = new uint64_t[size+1];
  uint8_t *reverseH = new uint8_t[size];
  size_t reverseOrderPos;

  // the lowest 2 bits are the h index (0, 1, or 2)
  // so we only have 6 bits for counting;
  // but that's sufficient
  uint8_t *t2count = new uint8_t[arrayLength];
  
  uint64_t *t2hash = new uint64_t[arrayLength];

  size_t *alone = new size_t[arrayLength];
  hashIndex = 0;
  
  // the array h0, h1, h2, h0, h1, h2
  size_t h012[5];

  while (true) {
    memset(t2count, 0, sizeof(uint8_t[arrayLength]));
    memset(t2hash, 0, sizeof(uint64_t[arrayLength]));

    // counting sort

    memset(reverseOrder, 0, sizeof(uint64_t[size]));
    reverseOrder[size] = 1;

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
      while(reverseOrder[startPos[segment_index]] != 0) {
        segment_index++;
        segment_index &= (size_t(1) << blockBits) - 1;
      }
      reverseOrder[startPos[segment_index]] = hash;
      startPos[segment_index]++;
    }
    uint8_t countMask = 0;
    for (size_t i = 0; i < size; i++) {
      uint64_t hash = reverseOrder[i];
      for (int hi = 0; hi < 3; hi++) {
        int index = getHashFromHash(hash, hi);
        t2count[index] += 4;
        t2count[index] ^= hi;
        t2hash[index] ^= hash;
        countMask |= t2count[index];
      }
    }
    delete[] startPos;
    if (countMask >= 0x80) {
      // we have a possible counter overflow
      // this branch is never taken except if there is a problem in the hash code
      // in which case construction fails
      continue;
    }

    reverseOrderPos = 0;
    size_t alonePos = 0;
    for (size_t i = 0; i < arrayLength; i++) {
      alone[alonePos] = i;
      int inc = (t2count[i] >> 2) == 1 ? 1 : 0;
      alonePos += inc;
    }

    while (alonePos > 0) {
      alonePos--;
      size_t index = alone[alonePos];
      if ((t2count[index] >> 2) == 1) {
        // It is still there!
        uint64_t hash = t2hash[index];
        int found = t2count[index] & 3;
        
        reverseH[reverseOrderPos] = found;
        reverseOrder[reverseOrderPos] = hash;
        
        h012[0] = getHashFromHash(hash, 0);
        h012[1] = getHashFromHash(hash, 1);
        h012[2] = getHashFromHash(hash, 2);

        size_t index3 = h012[mod3(found + 1)];
        alone[alonePos] = index3;
        alonePos += ((t2count[index3] >> 2) == 2 ? 1 : 0);
        t2count[index3] -= 4;
        t2count[index3] ^= mod3(found + 1);
        t2hash[index3] ^= hash;

        index3 = h012[mod3(found + 2)];
        alone[alonePos] = index3;
        alonePos += ((t2count[index3] >> 2) == 2 ? 1 : 0);
        t2count[index3] -= 4;
        t2count[index3] ^= mod3(found + 2);
        t2hash[index3] ^= hash;

        reverseOrderPos++;
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
  delete[] t2count;
  delete[] t2hash;
  
  for (int i = reverseOrderPos - 1; i >= 0; i--) {
    uint64_t hash = reverseOrder[i];
    int found = reverseH[i];
    FingerprintType xor2 = fingerprint(hash);
    h012[0] = getHashFromHash(hash, 0);
    h012[1] = getHashFromHash(hash, 1);
    h012[2] = getHashFromHash(hash, 2);
    h012[3] = h012[0];
    h012[4] = h012[1];
    fingerprints[h012[found]] = xor2 ^ fingerprints[h012[found + 1]] ^ fingerprints[h012[found + 2]];
  }
  delete[] reverseOrder;
  delete[] reverseH;

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
} // namespace xorbinaryfusefilter_lowmem
#endif