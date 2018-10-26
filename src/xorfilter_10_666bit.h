#ifndef XOR_FILTER10_666_XOR_FILTER_H_
#define XOR_FILTER10_666_XOR_FILTER_H_

#include <assert.h>
#include <algorithm>
#include <xorfilter.h>
#include <xorfilter_10bit.h>
#include "hashutil.h"

using namespace std;
using namespace hashing;

namespace xorfilter {

const uint32_t fingerMul = 1625;
// const uint32_t fingerMul = 1024;
const uint32_t fingerMul2 = fingerMul * fingerMul;
// (1<<64) / fingerMul
const __uint128_t invFingerMul = 11351842506898185L;
// (1<<64) / fingerMul2
const __uint128_t invFingerMul2 = 6985749235014L;

template <typename ItemType, typename HashFamily = TwoIndependentMultiplyShift>
class XorFilter10_666 {

  size_t size;
  size_t arrayLength;
  size_t blockLength;
  uint32_t *fingerprints;

  HashFamily* hasher;

  inline uint32_t fingerprint(const uint64_t hash) const {
    return (uint32_t) reduce(hash ^ (hash >> 32), fingerMul);
  }

 public:
  explicit XorFilter10_666(const size_t size) {
    hasher = new HashFamily();
    this->size = size;
    this->arrayLength = 32 + 1.23 * size;
    this->blockLength = arrayLength / 3;
    fingerprints = new uint32_t[blockLength]();
  }

  ~XorFilter10_666() {
    delete[] fingerprints;
    delete hasher;
  }

  Status AddAll(const vector<ItemType> data, const size_t start, const size_t end);

  // Report if the item is inserted, with false positive rate.
  Status Contain(const ItemType &item) const;

  // number of current inserted items;
  size_t Size() const { return size; }

  // size of the filter in bytes.
  size_t SizeInBytes() const { return blockLength * sizeof(uint32_t); }
};

int getFinger(uint32_t x, int i) {
    if (i == 0) {
        return x;
    }
    if (i == 1) {
//        return x >> 10;
        return (int) (((__uint128_t) x * (invFingerMul + 1)) >> 64);
    }
    if (i == 2) {
//        return x >> 20;
        return (int) (((__uint128_t) x * (invFingerMul2 + 1)) >> 64);
    }
    exit(0);
}

template <typename ItemType, typename HashFamily>
Status XorFilter10_666<ItemType, HashFamily>::AddAll(
    const vector<ItemType> keys, const size_t start, const size_t end) {

    int m = arrayLength;
    uint64_t* reverseOrder = new uint64_t[size];
    uint8_t* reverseH = new uint8_t[size];
    size_t reverseOrderPos;
    int hashIndex = 0;
    t2val_t * t2vals = new t2val_t[m];
    while (true) {
        memset(t2vals, 0, sizeof(t2val_t[m]));
        int blocks = 1 + ((3 * blockLength) >> blockShift);
        uint64_t* tmp = new uint64_t[blocks << blockShift];
        int* tmpc = new int[blocks]();
        for(size_t i = start; i < end; i++) {
            uint64_t k = keys[i];
            uint64_t hash = (*hasher)(k);
            for (int hi = 0; hi < 3; hi++) {
                int index = getHashFromHash(hash, hi, blockLength);
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
                int h = getHashFromHash(hash, hi, blockLength);
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

/*

        int* alone = new int[arrayLength];
        int alonePos = 0;
        for (size_t i = 0; i < arrayLength; i++) {
            if (t2vals[i].t2count == 1) {
                alone[alonePos++] = i;
            }
        }
        reverseOrderPos = 0;
        while (alonePos > 0 && reverseOrderPos < size) {
            int i = alone[--alonePos];
            if (t2vals[i].t2count == 0) {
                continue;
            }
            long hash = t2vals[i].t2;
            uint8_t found = -1;
            for (int hi = 0; hi < 3; hi++) {
                int h = getHashFromHash(hash, hi, blockLength);
                int newCount =  --t2vals[h].t2count;
                if (newCount == 0) {
                    found = (uint8_t) hi;
                } else {
                    if (newCount == 1) {
                        alone[alonePos++] = h;
                    }
                    t2vals[h].t2 ^= hash;
                }
            }
            reverseOrder[reverseOrderPos] = hash;
            reverseH[reverseOrderPos] = found;
            reverseOrderPos++;
        }
        delete [] alone;
*/

        if (reverseOrderPos == size) {
            break;
        }

        std::cout << "WARNING: hashIndex " << hashIndex << "\n";
        if (hashIndex >= 0) {
            std::cout << (end - start) << " keys; arrayLength " << arrayLength
                << " blockLength " << blockLength
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
        uint32_t sum = fingerprint(hash);
        for (int hi = 0; hi < 3; hi++) {
            size_t h = getHashFromHash10(hash, hi, blockLength);
            if (found == hi) {
                change = h;
            } else {
                // this is different from BDZ: using xor to calculate the
                // fingerprint
                sum += getFinger(fingerprints[h], hi) % fingerMul;
            }
        }
        // set such that if the right fingerprint is added, the result is 0
        uint32_t set = (fingerMul - (sum % fingerMul)) % fingerMul;
        fingerprints[change] += found == 0 ? set : found == 1 ? (set * fingerMul) : (set * fingerMul2);
    }
    delete [] t2vals;
    delete [] reverseOrder;
    delete [] reverseH;

    return Ok;
}

template <typename ItemType, typename HashFamily>
Status XorFilter10_666<ItemType, HashFamily>::Contain(
    const ItemType &key) const {
    uint64_t hash = (*hasher)(key);
    uint32_t f = fingerprint(hash);
    uint32_t r0 = (uint32_t) hash;
    uint32_t r1 = (uint32_t) rotl64(hash, 21);
    uint32_t r2 = (uint32_t) rotl64(hash, 42);
    uint32_t h0 = reduce(r0, blockLength);
    uint32_t h1 = reduce(r1, blockLength);
    uint32_t h2 = reduce(r2, blockLength);
    f += fingerprints[h0];
    f += (((__uint128_t) fingerprints[h1] * (invFingerMul + 1)) >> 64);
    f += (((__uint128_t) fingerprints[h2] * (invFingerMul2 + 1)) >> 64);
    return (f % fingerMul) == 0 ? Ok : NotFound;
}

}  // namespace xorfilter
#endif  // XOR_FILTER10_666_XOR_FILTER_H_
