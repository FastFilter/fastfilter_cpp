#ifndef XOR_FUSE_FILTER_XOR_FILTER_H_
#define XOR_FUSE_FILTER_XOR_FILTER_H_

#include <assert.h>
#include <algorithm>
#include "hashutil.h"

using namespace std;
using namespace hashing;

namespace xorfusefilter {
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
    return (n << c) | ( n >> ((-c) & mask));
}

__attribute__((always_inline))
inline uint32_t reduce(uint32_t hash, uint32_t n) {
    // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
    return (uint32_t) (((uint64_t) hash * n) >> 32);
}

#define H128
const size_t segmentLengthBits = 13;
const size_t segmentLength = 1 << segmentLengthBits;

size_t getHashFromHash(uint64_t hash, int index, int segmentCount) {
#ifdef H128
    __uint128_t x = (__uint128_t) hash * (__uint128_t) segmentCount;
    int seg = (uint64_t)(x >> 64);
    uint64_t hh = hash;
#else
    uint64_t seg = reduce(hash, segmentCount);
    uint64_t hh = (hash ^ (hash >> 32));
#endif
    int h = (seg + index) * segmentLength + (size_t)((hh >> (index * segmentLengthBits)) & (segmentLength - 1));
    return h;
}

template <typename ItemType, typename FingerprintType,
          typename HashFamily = TwoIndependentMultiplyShift>
class XorFuseFilter {
 public:

  size_t size;
  size_t arrayLength;
  size_t segmentCount;
  FingerprintType *fingerprints;

  HashFamily* hasher;

  inline FingerprintType fingerprint(const uint64_t hash) const {
    return (FingerprintType) hash;
    // return (FingerprintType) hash ^ (hash >> 32);
  }

  explicit XorFuseFilter(const size_t size) {
    hasher = new HashFamily();
    this->size = size;
    size_t capacity = size / 0.879;
    capacity = (capacity + 3 - 1) / 3 * 3;
    capacity = (capacity + segmentLength - 1) / segmentLength * segmentLength;
    this->segmentCount = capacity / segmentLength;
    this->arrayLength = (segmentCount + 2) * segmentLength;
    fingerprints = new FingerprintType[arrayLength]();
    std::fill_n(fingerprints, arrayLength, 0);
  }

  ~XorFuseFilter() {
    delete[] fingerprints;
    delete hasher;
  }

  Status AddAll(const vector<ItemType> &data, const size_t start, const size_t end) {
      return AddAll(data.data(),start,end);
  }

  Status AddAll(const ItemType* data, const size_t start, const size_t end);

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

const int blockShift = 18;

void applyBlock(uint64_t* tmp, int b, int len, t2val_t * t2vals) {
    for (int i = 0; i < len; i += 2) {
        uint64_t x = tmp[(b << blockShift) + i];
        int index = (int) tmp[(b << blockShift) + i + 1];
        t2vals[index].t2count++;
        t2vals[index].t2 ^= x;
    }
}

int applyBlock2(uint64_t* tmp, int b, int len, t2val_t * t2vals, int* alone, int alonePos) {
    for (int i = 0; i < len; i += 2) {
        uint64_t hash = tmp[(b << blockShift) + i];
        int index = (int) tmp[(b << blockShift) + i + 1];
        int oldCount = t2vals[index].t2count;
        if (oldCount >= 1) {
            int newCount = oldCount - 1;
            t2vals[index].t2count = newCount;
            if (newCount == 1) {
                alone[alonePos++] = index;
            }
            t2vals[index].t2 ^= hash;
        }
    }
    return alonePos;
}

template <typename ItemType, typename FingerprintType,
          typename HashFamily>
Status XorFuseFilter<ItemType, FingerprintType, HashFamily>::AddAll(
    const ItemType* keys, const size_t start, const size_t end) {

    int m = arrayLength;
    uint64_t* reverseOrder = new uint64_t[size];
    uint8_t* reverseH = new uint8_t[size];
    size_t reverseOrderPos;
    int hashIndex = 0;
    t2val_t * t2vals = new t2val_t[m];
    while (true) {
        memset(t2vals, 0, sizeof(t2val_t[m]));
        int blocks = 1 + (arrayLength >> blockShift);
        uint64_t* tmp = new uint64_t[blocks << blockShift];
        int* tmpc = new int[blocks]();
        for(size_t i = start; i < end; i++) {
            uint64_t k = keys[i];
            uint64_t hash = (*hasher)(k);
            for (int hi = 0; hi < 3; hi++) {
                int index = getHashFromHash(hash, hi, segmentCount);
                int b = index >> blockShift;
                int i2 = tmpc[b];
                tmp[(b << blockShift) + i2] = hash;
                tmp[(b << blockShift) + i2 + 1] = index;
                tmpc[b] += 2;
                if (i2 + 2 == (1 << blockShift)) {
                    applyBlock(tmp, b, i2 + 2, t2vals);
                    tmpc[b] = 0;
                }
            }

        }
        for (int b = 0; b < blocks; b++) {
            applyBlock(tmp, b, tmpc[b], t2vals);
        }
        delete[] tmp;
        delete[] tmpc;
        reverseOrderPos = 0;

        int* alone = new int[arrayLength];
        int alonePos = 0;
        for (size_t i = 0; i < arrayLength; i++) {
            if (t2vals[i].t2count == 1) {
                alone[alonePos++] = i;
            }
        }
        tmp = new uint64_t[blocks << blockShift];
        tmpc = new int[blocks]();
        reverseOrderPos = 0;
        int bestBlock = -1;
        while (reverseOrderPos < size) {
            if (alonePos == 0) {
                // we need to apply blocks until we have an entry that is alone
                // (that is, until alonePos > 0)
                // so, find a large block (the larger the better)
                // but don't need to search very long
                // start searching where we stopped the last time
                // (to make it more even)
                for (int i = 0, b = bestBlock + 1, best = -1; i < blocks; i++) {
                    if (b >= blocks) {
                        b = 0;
                    }
                    if (tmpc[b] > best) {
                        best = tmpc[b];
                        bestBlock = b;
                        if (best > (1 << (blockShift - 1))) {
                            // sufficiently large: stop
                            break;
                        }
                    }
                }
                if (tmpc[bestBlock] > 0) {
                    alonePos = applyBlock2(tmp, bestBlock, tmpc[bestBlock], t2vals, alone, alonePos);
                    tmpc[bestBlock] = 0;
                }
                // applying a block may not actually result in a new entry that is alone
                if (alonePos == 0) {
                    for (int b = 0; b < blocks && alonePos == 0; b++) {
                        if (tmpc[b] > 0) {
                            alonePos = applyBlock2(tmp, b, tmpc[b], t2vals, alone, alonePos);
                            tmpc[b] = 0;
                        }
                    }
                }
            }
            if (alonePos == 0) {
                break;
            }
            int i = alone[--alonePos];
            int b = i >> blockShift;
            if (tmpc[b] > 0) {
                alonePos = applyBlock2(tmp, b, tmpc[b], t2vals, alone, alonePos);
                tmpc[b] = 0;
            }
            uint8_t found = -1;
            if (t2vals[i].t2count == 0) {
                continue;
            }
            long hash = t2vals[i].t2;
            for (int hi = 0; hi < 3; hi++) {
                int h = getHashFromHash(hash, hi, segmentCount);
                if (h == i) {
                    found = (uint8_t) hi;
                    t2vals[i].t2count = 0;
                } else {
                    int b = h >> blockShift;
                    int i2 = tmpc[b];
                    tmp[(b << blockShift) + i2] = hash;
                    tmp[(b << blockShift) + i2 + 1] = h;
                    tmpc[b] += 2;
                    if (tmpc[b] >= 1 << blockShift) {
                        alonePos = applyBlock2(tmp, b, tmpc[b], t2vals, alone, alonePos);
                        tmpc[b] = 0;
                    }
                }
            }
            reverseOrder[reverseOrderPos] = hash;
            reverseH[reverseOrderPos] = found;
            reverseOrderPos++;
        }
        delete[] tmp;
        delete[] tmpc;
        delete[] alone;

        if (reverseOrderPos == size) {
            break;
        }

        std::cout << "WARNING: hashIndex " << hashIndex << "\n";
        if (hashIndex >= 0) {
            std::cout << (end - start) << " keys; arrayLength " << arrayLength
                << " reverseOrderPos " << reverseOrderPos << "\n";
        }

        hashIndex++;

        // use a new random numbers
        delete hasher;
        hasher = new HashFamily();

    }

    for (int i = reverseOrderPos - 1; i >= 0; i--) {
        // the hash of the key we insert next
        uint64_t hash = reverseOrder[i];
        int found = reverseH[i];
        // which entry in the table we can change
        int change = -1;
        // we set table[change] to the fingerprint of the key,
        // unless the other two entries are already occupied
        FingerprintType xor2 = fingerprint(hash);
        for (int hi = 0; hi < 3; hi++) {
            size_t h = getHashFromHash(hash, hi, segmentCount);
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
    delete [] t2vals;
    delete [] reverseOrder;
    delete [] reverseH;

    return Ok;
}

template <typename ItemType, typename FingerprintType,
          typename HashFamily>
Status XorFuseFilter<ItemType, FingerprintType, HashFamily>::Contain(
    const ItemType &key) const {
    uint64_t hash = (*hasher)(key);
    FingerprintType f = fingerprint(hash);
#ifdef H128
    __uint128_t x = (__uint128_t) hash * (__uint128_t) segmentCount;
    int seg = (uint64_t)(x >> 64);
    uint64_t hh = hash;
#else
    uint64_t seg = reduce(hash, segmentCount);
    uint64_t hh = (hash ^ (hash >> 32));
#endif
    int h0 = (seg + 0) * segmentLength + (size_t)((hh >> (0 * segmentLengthBits)) & (segmentLength - 1));
    int h1 = (seg + 1) * segmentLength + (size_t)((hh >> (1 * segmentLengthBits)) & (segmentLength - 1));
    int h2 = (seg + 2) * segmentLength + (size_t)((hh >> (2 * segmentLengthBits)) & (segmentLength - 1));
    f ^= fingerprints[h0] ^ fingerprints[h1] ^ fingerprints[h2];
    return f == 0 ? Ok : NotFound;
}

template <typename ItemType, typename FingerprintType,
          typename HashFamily>
std::string XorFuseFilter<ItemType, FingerprintType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "XorFuseFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  return ss.str();
}
}  // namespace xorfusefilter
#endif  // XOR_FUSE_FILTER_XOR_FILTER_H_
