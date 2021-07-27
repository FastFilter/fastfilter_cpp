#ifndef FOURWISE_XOR_BINARY_FUSE_FILTER_XOR_FILTER_NAIVE_H_
#define FOURWISE_XOR_BINARY_FUSE_FILTER_XOR_FILTER_NAIVE_H_
#include "xor_binary_fuse_filter.h"
namespace xorbinaryfusefilter_naive4wise {
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
  static constexpr size_t arity = 4;
  FingerprintType *fingerprints;

  HashFamily *hasher;

  inline FingerprintType fingerprint(const uint64_t hash) const {
    return (FingerprintType)hash;
  }

  static inline __attribute__((always_inline)) uint64_t rotate(uint64_t n,
                                                               unsigned int c) {
    const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
    c &= mask;
    return (n << c) | (n >> ((-c) & mask));
  }

  inline __attribute__((always_inline)) size_t
  getHashFromHash(uint64_t hash, int index) const {
    __uint128_t x = (__uint128_t)hash * (__uint128_t)segmentCountLength;
    uint64_t h = (uint64_t)(x >> 64);
    h += index * segmentLength;
    uint64_t hh = hash;
    h ^= rotate(hh, index * 16) & segmentLengthMask;
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
  uint64_t t2;
  uint64_t t2count;
};

typedef struct t2val t2val_t;

template <typename ItemType, typename FingerprintType, typename HashFamily>
Status XorBinaryFuseFilter<ItemType, FingerprintType, HashFamily>::AddAll(
    const ItemType *keys, const size_t start, const size_t end) {

  uint64_t *reverseOrder = new uint64_t[size];
  uint8_t *reverseH = new uint8_t[size];
  size_t reverseOrderPos;

  t2val_t *t2vals = new t2val_t[arrayLength];

  size_t *alone = new size_t[arrayLength];
  size_t hashIndex{0};

  while (true) {
    memset(t2vals, 0, sizeof(t2val_t[arrayLength]));
    for (size_t i = start; i < end; i++) {
      uint64_t k = keys[i];
      uint64_t hash = (*hasher)(k);
      for (int hi = 0; hi < 4; hi++) {
        int index = getHashFromHash(hash, hi);
        t2vals[index].t2count++;
        t2vals[index].t2 ^= hash;
      }
    }
    reverseOrderPos = 0;
    size_t alonePos = 0;
    for (size_t i = 0; i < arrayLength; i++) {
      if (t2vals[i].t2count == 1) {
        alone[alonePos++] = i;
      }
    }
    while (alonePos > 0) {
      alonePos--;
      size_t index = alone[alonePos];
      if (t2vals[index].t2count == 1) {
        // It is still there!
        uint64_t hash = t2vals[index].t2;
        reverseOrder[reverseOrderPos] = hash;
        for (int hi = 0; hi < 4; hi++) {
          size_t index3 = getHashFromHash(hash, hi);
          if (index3 == index) {
            reverseH[reverseOrderPos] = hi;
            // no need to decrement & remove
            continue;
          } else if (t2vals[index3].t2count == 2) {
            // Found a new candidate !
            alone[alonePos++] = index3;
          }
          t2vals[index3].t2count -= 1;
          t2vals[index3].t2 ^= hash;
        }
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
  delete[] t2vals;

  for (int i = reverseOrderPos - 1; i >= 0; i--) {
    // the hash of the key we insert next
    uint64_t hash = reverseOrder[i];
    int found = reverseH[i];
    // which entry in the table we can change
    int change = -1;
    // we set table[change] to the fingerprint of the key,
    // unless the other two entries are already occupied
    FingerprintType xor2 = fingerprint(hash);
    for (int hi = 0; hi < 4; hi++) {
      size_t h = getHashFromHash(hash, hi);
      if (found == hi) {
        change = h;
      } else {
        // this is different from BDZ: using xor to calculate the
        // fingerprint
        xor2 ^= fingerprints[h];
      }
    }
    fingerprints[change] = xor2;
  }
  delete[] reverseOrder;
  delete[] reverseH;

  return Ok;
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
Status XorBinaryFuseFilter<ItemType, FingerprintType, HashFamily>::Contain(
    const ItemType &key) const {
  uint64_t hash = (*hasher)(key);
  // Could manually optimize.
  FingerprintType f = fingerprint(hash);
  for (int hi = 0; hi < 4; hi++) {
    size_t h = getHashFromHash(hash, hi);
    f ^= fingerprints[h];
  }
  return f == 0 ? Ok : NotFound;
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
std::string
XorBinaryFuseFilter<ItemType, FingerprintType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "4-wise XorBinaryFuseFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  return ss.str();
}
} // namespace xorbinaryfusefilter_naive4wise

#endif