#ifndef XOR_FILTER2N_XOR_FILTER2N_H_
#define XOR_FILTER2N_XOR_FILTER2N_H_

#include <assert.h>
#include <algorithm>

#include "debug.h"
#include "hashutil.h"
#include "printutil.h"
#include "nbit_array.h"

using namespace std;
using namespace cuckoofilter;

namespace xorfilter2n {
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

size_t getHashFromHash(uint64_t hash, int index, int blockLength) {
    uint32_t r;
    switch(index) {
    case 0:
        r = (uint32_t) (hash);
        break;
    case 1:
        r = (uint32_t) rotl64(hash, 21);
        break;
    default:
        r = (uint32_t) rotl64(hash, 42);
        break;
    }
    r = r & (blockLength - 1);
    r = r + index * blockLength;
    return (size_t) r;
}

template <typename ItemType, typename FingerprintType,
          typename FingerprintStorageType, typename HashFamily = TwoIndependentMultiplyShift>
class XorFilter2n {

  size_t size;
  size_t arrayLength;
  size_t blockLength;
  size_t blockMask;
  FingerprintStorageType *fingerprints;
  uint64_t fingerprintMask;

  HashFamily* hasher;

  inline FingerprintType fingerprint(const uint64_t hash) const {
    return (FingerprintType) fingerprints->mask(hash);
  }

 public:
  explicit XorFilter2n(const size_t size) {
    hasher = new HashFamily();
    this->size = size;
    this->arrayLength = 3 + 1.23 * size;
    this->blockLength = 1;
    while (this->blockLength < arrayLength / 3) {
        this->blockLength *= 2;
    }
    blockMask = blockLength - 1;
    this->arrayLength = 3 * blockLength;
    fingerprints = new FingerprintStorageType(arrayLength);
  }

  ~XorFilter2n() {
    delete fingerprints;
    delete hasher;
  }

  Status AddAll(const vector<ItemType> data, const size_t start, const size_t end);

  // Report if the item is inserted, with false positive rate.
  Status Contain(const ItemType &item) const;

  /* methods for providing stats  */
  // summary infomation
  std::string Info() const;

  // number of current inserted items;
  size_t Size() const { return size; }

  // size of the filter in bytes.
  size_t SizeInBytes() const { return fingerprints->getByteCount(); }
};

template <typename ItemType, typename FingerprintType,
          typename FingerprintStorageType, typename HashFamily>
Status XorFilter2n<ItemType, FingerprintType, FingerprintStorageType, HashFamily>::AddAll(
    const vector<ItemType> keys, const size_t start, const size_t end) {
    int m = arrayLength;
    uint64_t* reverseOrder = new uint64_t[size];
    uint8_t* reverseH = new uint8_t[size];
    size_t reverseOrderPos;
    int hashIndex = 0;
    uint64_t* t2 = new uint64_t[m];
    uint8_t* t2count = new uint8_t[m];
    while (true) {
        memset(t2count, 0, sizeof(uint8_t[m]));
        memset(t2, 0, sizeof(uint64_t[m]));
        for(size_t i = start; i < end; i++) {
            uint64_t k = keys[i];
            uint64_t hash = (*hasher)(k);
            int h0 = (int) (hash & blockMask);
            int h1 = (int) (rotl64(hash, 21) & blockMask) + blockLength;
            int h2 = (int) (rotl64(hash, 42) & blockMask) + 2 * blockLength;
            t2count[h0]++;
            t2[h0] ^= hash;
            t2count[h1]++;
            t2[h1] ^= hash;
            t2count[h2]++;
            t2[h2] ^= hash;
        }
        reverseOrderPos = 0;
        int* alone = new int[arrayLength];
        int alonePos = 0;
        reverseOrderPos = 0;
        for(size_t nextAloneCheck = 0; nextAloneCheck < arrayLength;) {
            while (nextAloneCheck < arrayLength) {
                if (t2count[nextAloneCheck] == 1) {
                    alone[alonePos++] = nextAloneCheck;
                    // break;
                }
                nextAloneCheck++;
            }
            while (alonePos > 0) {
                int i = alone[--alonePos];
                if (t2count[i] == 0) {
                    continue;
                }
                long hash = t2[i];
                uint8_t found = -1;
                for (int hi = 0; hi < 3; hi++) {
                    int h = getHashFromHash(hash, hi, blockLength);
                    int newCount = --t2count[h];
                    if (newCount == 0) {
                        found = (uint8_t) hi;
                    } else {
                        if (newCount == 1) {
                            alone[alonePos++] = h;
                        }
                        t2[h] ^= hash;
                    }
                }
                reverseOrder[reverseOrderPos] = hash;
                reverseH[reverseOrderPos] = found;
                reverseOrderPos++;
            }
        }
        delete [] alone;
        if (reverseOrderPos == size) {
            break;
        }

        std::cout << "WARNING: hashIndex " << hashIndex << "\n";
        if (hashIndex >= 0) {
            size_t outputlimit = 5; // we don't want to spam
            std::cout << (end - start) << " keys; arrayLength " << arrayLength
                << " blockLength " << blockLength
                << " reverseOrderPos " << reverseOrderPos << "\n";

            int pos = 0;
            for (size_t i = 0; pos < 1000 && i < arrayLength; i++) {
                if (t2count[i] > 1) {
                    if(outputlimit > 0) {
                       std::cout << "  count[" << i << "] = " << (int) t2count[i] << "\n";
                       outputlimit --;
                     }
                }
            }
           for(size_t i = start; i < end; i++) {
                uint64_t k = keys[i];
                uint64_t hash = (*hasher)(k);
                int h0 = (int) (hash & blockMask);
                int h1 = (int) (rotl64(hash, 21) & blockMask) + blockLength;
                int h2 = (int) (rotl64(hash, 42) & blockMask) + 2 * blockLength;
                if (t2count[h0] > 1 || t2count[h1] > 1 || t2count[h2] > 1) {
                    if(outputlimit > 0) {
                      std::cout << "  key " << k << " hash=" << hash << " h0=" << h0 << " h1=" << h1 << " h2=" << h2 << "\n";
                      outputlimit --;
                    }
                }
            }

            // for(size_t i = start; i < end; i++) {
            //     uint64_t k = keys[i];
            //     std::cout << k << "\n";
            // }
            // std::cout << "end\n";
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
            size_t h = getHashFromHash(hash, hi, blockLength);
            if (found == hi) {
                change = h;
            } else {
                // this is different from BDZ: using xor to calculate the
                // fingerprint
                xor2 ^= fingerprints->mask(fingerprints->get(h));
            }
        }
        fingerprints->set(change, xor2);
    }

    delete [] t2;
    delete [] t2count;
    delete [] reverseOrder;
    delete [] reverseH;

    return Ok;
}

template <typename ItemType, typename FingerprintType,
          typename FingerprintStorageType, typename HashFamily>
Status XorFilter2n<ItemType, FingerprintType, FingerprintStorageType, HashFamily>::Contain(
    const ItemType &key) const {
    uint64_t hash = (*hasher)(key);
    FingerprintType f = hash;
    uint32_t r0 = (uint32_t) hash;
    uint32_t r1 = (uint32_t) rotl64(hash, 21);
    uint32_t r2 = (uint32_t) rotl64(hash, 42);
    uint32_t h0 = (r0 & blockMask);
    uint32_t h1 = (r1 & blockMask) + blockLength;
    uint32_t h2 = (r2 & blockMask) + 2 * blockLength;
    f ^= fingerprints->get(h0) ^ fingerprints->get(h1) ^ fingerprints->get(h2);
    return fingerprint(f) == 0 ? Ok : NotFound;
}

template <typename ItemType, typename FingerprintType,
          typename FingerprintStorageType, typename HashFamily>
std::string XorFilter2n<ItemType, FingerprintType, FingerprintStorageType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "XorFilter2n Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  return ss.str();
}
}  // namespace xorfilter2n
#endif  // XOR_FILTER2N_XOR_FILTER2N_H_
