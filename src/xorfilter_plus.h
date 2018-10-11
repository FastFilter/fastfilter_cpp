#ifndef XOR_FILTER_PLUS_XOR_FILTER_PLUS_H_
#define XOR_FILTER_PLUS_XOR_FILTER_PLUS_H_

#include <assert.h>
#include <algorithm>

#include "debug.h"
#include "hashutil.h"
#include "printutil.h"

using namespace std;
using namespace cuckoofilter;

namespace xorfilter_plus {
// status returned by a xor filter operation
enum Status {
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

inline int numberOfLeadingZeros64(uint64_t x) {
    // If x is 0, the result is undefined.
    return __builtin_clzl(x);
}

inline int mostSignificantBit(uint64_t x) {
    return 63 - numberOfLeadingZeros64(x);
}

inline int bitCount64(uint64_t x) {
    return __builtin_popcountll(x);
}

class Rank9 {

    uint64_t* bits;
    uint64_t bitsArraySize;
    uint64_t* counts;
    uint64_t countsArraySize;

public:

    Rank9(uint64_t* sourceBits, size_t bitCount) {
        // One zero entry is needed at the end
        bitsArraySize = 1 + (size_t) ((bitCount + 63) / 64);
        bits = new uint64_t[bitsArraySize];
        memcpy(bits, sourceBits, (bitsArraySize - 1) * sizeof(uint64_t));
        bits[bitsArraySize - 1] = 0;
        uint64_t length = bitsArraySize * 64;
        size_t numWords = (size_t) ((length + 63) / 64);
        size_t numCounts = (size_t) ((length + 8 * 64 - 1) / (8 * 64)) * 2;
        countsArraySize = numCounts + 1;
        counts = new uint64_t[countsArraySize];
        // just to be sure
        memset(counts, 0, sizeof(uint64_t[countsArraySize]));
        uint64_t c = 0;
        uint64_t pos = 0;
        for (uint64_t i = 0; i < numWords; i += 8, pos += 2) {
            counts[pos] = c;
            counts[pos + 1] = 0;
            c += bitCount64(bits[i]);
            for (uint64_t j = 1; j < 8; j++) {
                counts[pos + 1] |= (c - counts[pos]) << 9 * (j - 1);
                if (i + j < numWords) {
                    c += bitCount64(bits[i + j]);
                }
            }
        }
        counts[numCounts] = c;
        /*
        // self-test (should be somewhere else)
        for(int i=0, j = 0; i<bitCount;i++) {
            uint64_t bit1 = sourceBits[i >> 6] & (1L << (i & 63));
            uint64_t bit2 = (bits[i >> 6] >> (i & 63)) & 1;
            if ((bit2 == 0 ? 0 : 1) != (bit1 == 0 ? 0 : 1)) {
                std::cout << "WARNING: get incorrect " << i << "\n";
            }
            if (rank(i) != j) {
                std::cout << "WARNING: rank incorrect " << i << "\n";
            }
            int rag = getAndPartialRank(i) + (remainingRank(i) << 1);
            if ((rag >> 1) != j) {
                std::cout << "WARNING: getAndPartialRank incorrect " << i << "\n";
            }
            if (bit2 != 0) {
                j++;
            }
        }
        */
    }

    ~Rank9() {
        delete[] bits;
        delete[] counts;
    }

    uint64_t rank(uint64_t pos) {
        uint64_t word = pos >> 6;
        uint64_t block = (word >> 2) & ~1;
        int32_t offset = (word & 7) - 1;
        return counts[block] +
                ((counts[block + 1] >> (offset + ((offset >> 28) & 8)) * 9) & 0x1ff) +
                bitCount64(bits[word] & ((1L << (pos & 63)) - 1));
    }

    uint64_t get(uint64_t pos) {
        return (bits[(size_t) (pos >> 6)] >> pos) & 1;
    }

    uint64_t getAndPartialRank(uint64_t pos) {
        uint64_t word = pos >> 6;
        uint64_t x = bits[word];
        return ((bitCount64(x & ((1L << (pos & 63)) - 1))) << 1) +
                ((x >> (pos & 63)) & 1);
    }

    uint64_t remainingRank(uint64_t pos) {
        uint64_t word = pos >> 6;
        uint64_t block = (word >> 2) & ~1;
        int32_t offset = (word & 7) - 1;
        return counts[block] +
                ((counts[block + 1] >> (offset + ((offset >> 28) & 8)) * 9) & 0x1ff);
    }

    uint64_t getBitCount() {
        return bitsArraySize * 64 + countsArraySize * 64;
    }

};

inline uint64_t rotl64(uint64_t n, unsigned int c) {
    // assumes width is a power of 2
    const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
    // assert ( (c<=mask) &&"rotate by type width or more");
    c &= mask;
    return (n << c) | ( n >> ((-c) & mask));
}

inline uint32_t reduce(uint32_t hash, uint32_t n) {
    // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
    return (uint32_t) (((uint64_t) hash * n) >> 32);
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
    r = reduce(r, blockLength);
    r = r + index * blockLength;
    return (size_t) r;
}

struct t2val {
  uint64_t t2;
  uint64_t t2count;
};

typedef struct t2val t2val_t;

#define BLOCK_SHIFT 18
#define BLOCK_LEN (1 << BLOCK_SHIFT)

void applyBlock(uint64_t* tmp, int b, int len, t2val_t * t2vals) {
    for (int i = 0; i < len; i += 2) {
        uint64_t x = tmp[(b << BLOCK_SHIFT) + i];
        int index = (int) tmp[(b << BLOCK_SHIFT) + i + 1];
        t2vals[index].t2count++;
        t2vals[index].t2 ^= x;
    }
}

template <typename ItemType, typename FingerprintType,
          typename HashFamily = TwoIndependentMultiplyShift>
class XorFilterPlus {

  size_t size;
  size_t arrayLength;
  size_t blockLength;
  FingerprintType *fingerprints = NULL;
  Rank9 *rank = NULL;
  size_t totalSizeInBytes;

  HashFamily* hasher;

  inline FingerprintType fingerprint(const uint64_t hash) const {
    return (FingerprintType) hash;
  }

 public:
  explicit XorFilterPlus(const size_t size) {
    hasher = new HashFamily();
    this->size = size;
    this->arrayLength = 3 + 1.23 * size;
    this->blockLength = arrayLength / 3;
  }

  ~XorFilterPlus() {
    delete hasher;
    if (fingerprints != NULL) {
        delete[] fingerprints;
    }
    if (rank != 0) {
        delete rank;
    }
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
  size_t SizeInBytes() const { return totalSizeInBytes; }
};

template <typename ItemType, typename FingerprintType,
          typename HashFamily>
Status XorFilterPlus<ItemType, FingerprintType, HashFamily>::AddAll(
    const vector<ItemType> keys, const size_t start, const size_t end) {
    int m = arrayLength;
    uint64_t* reverseOrder = new uint64_t[size];
    uint8_t* reverseH = new uint8_t[size];
    size_t reverseOrderPos;
    int hashIndex = 0;
    t2val_t * t2vals = new t2val_t[m];
    while (true) {
        memset(t2vals, 0, sizeof(t2val_t[m]));
        int blocks = 1 + (3 * blockLength) / BLOCK_LEN;
        uint64_t* tmp = new uint64_t[blocks * BLOCK_LEN];
        int* tmpc = new int[blocks]();
        for(size_t i = start; i < end; i++) {
            uint64_t k = keys[i];
            uint64_t hash = (*hasher)(k);
            for (int hi = 0; hi < 3; hi++) {
                int index = getHashFromHash(hash, hi, blockLength);
                int b = index >> BLOCK_SHIFT;
                int i2 = tmpc[b];
                tmp[(b << BLOCK_SHIFT) + i2] = hash;
                tmp[(b << BLOCK_SHIFT) + i2 + 1] = index;
                tmpc[b] += 2;
                if (i2 + 2 == BLOCK_LEN) {
                    applyBlock(tmp, b, i2 + 2, t2vals);
                    tmpc[b] = 0;
                }
            }
/*
            int h0 = reduce((int) (hash), blockLength);
            int h1 = reduce((int) rotl64(hash, 21), blockLength) + blockLength;
            int h2 = reduce((int) rotl64(hash, 42), blockLength) + 2 * blockLength;
            t2vals[h0].t2count++;
            t2vals[h0].t2 ^= hash;
            t2vals[h1].t2count++;
            t2vals[h1].t2 ^= hash;
            t2vals[h2].t2count++;
            t2vals[h2].t2 ^= hash;
*/
        }
        for (int b = 0; b < blocks; b++) {
            applyBlock(tmp, b, tmpc[b], t2vals);
        }
        delete[] tmp;
        delete[] tmpc;

        reverseOrderPos = 0;
        int* alone[3];
        alone[0] = new int[blockLength];
        alone[1] = new int[blockLength];
        alone[2] = new int[blockLength];
        int alonePos[] = {0, 0, 0};
        for(int nextAlone = 0; nextAlone < 3; nextAlone++) {
            for (size_t i = 0; i < blockLength; i++) {
                if (t2vals[nextAlone * blockLength + i].t2count == 1) {
                    alone[nextAlone][alonePos[nextAlone]++] = nextAlone * blockLength + i;
                }
            }
        }
        int found = -1;
        while (true) {
            int i = -1;
            for (int hi = 0; hi < 3; hi++) {
                if (alonePos[hi] > 0) {
                    i = alone[hi][--alonePos[hi]];
                    found = hi;
                    break;
                }
            }
            if (i == -1) {
                // no entry found
                break;
            }
            if (t2vals[i].t2count <= 0) {
                continue;
            }
            uint64_t hash = t2vals[i].t2;
            --t2vals[i].t2count;
            // which index (0, 1, 2) the entry was found
            for (int hi = 0; hi < 3; hi++) {
                if (hi != found) {
                    int h = getHashFromHash(hash, hi, blockLength);
                    int newCount = --t2vals[h].t2count;
                    if (newCount == 1) {
                        // we found a key that is _now_ alone
                        alone[hi][alonePos[hi]++] = h;
                    }
                    // remove this key from the t2 table, using xor
                    t2vals[h].t2 ^= hash;
                }
            }
            reverseOrder[reverseOrderPos] = hash;
            reverseH[reverseOrderPos] = found;
            reverseOrderPos++;
        }
        delete [] alone[0];
        delete [] alone[1];
        delete [] alone[2];
        if (reverseOrderPos == size) {
            break;
        }

        std::cout << "WARNING: hashIndex " << hashIndex << "\n";
        if (hashIndex >= 0) {
           // size_t outputlimit = 5; // we don't want to spam
            std::cout << (end - start) << " keys; arrayLength " << arrayLength
                << " blockLength " << blockLength
                << " reverseOrderPos " << reverseOrderPos << "\n";
           // int pos = 0;
           /* for (size_t i = 0; pos < 1000 && i < arrayLength; i++) {
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
                int h0 = reduce((int) (hash), blockLength);
                int h1 = reduce((int) rotl64(hash, 21), blockLength) + blockLength;
                int h2 = reduce((int) rotl64(hash, 42), blockLength) + 2 * blockLength;
                if (t2count[h0] > 1 || t2count[h1] > 1 || t2count[h2] > 1) {
                    if(outputlimit > 0) {
                      std::cout << "  key " << k << " hash=" << hash << " h0=" << h0 << " h1=" << h1 << " h2=" << h2 << "\n";
                      outputlimit --;
                    }
                }
            }*/

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

    FingerprintType *fp = new FingerprintType[3 * blockLength];
    std::fill_n(fp, 3 * blockLength, 0);
    for (int i = reverseOrderPos - 1; i >= 0; i--) {
        // the hash of the key we insert next
        uint64_t hash = reverseOrder[i];
        int found = reverseH[i];
        // which entry in the table we can change
        int change = -1;
        // we set table[change] to the fingerprint of the key,
        // unless the other two entries are already occupied
        FingerprintType xor2 = (FingerprintType) fingerprint(hash);
        for (int hi = 0; hi < 3; hi++) {
            size_t h = getHashFromHash(hash, hi, blockLength);
            if (found == hi) {
                change = h;
            } else {
                // this is different from BDZ: using xor to calculate the
                // fingerprint
                xor2 ^= fp[h];
            }
        }
        fp[change] = xor2;
    }

    delete [] t2vals;
    delete [] reverseOrder;
    delete [] reverseH;

    uint64_t bitCount = blockLength;
    size_t wordlen = (bitCount + 63) / 64;
    uint64_t *bits = new uint64_t[wordlen];//was /63 but it is suspicious
    memset(bits,0,wordlen*sizeof(uint64_t)); // seems wrong that there is no init
    int setBits = 0;
    for (size_t i = 0; i < blockLength; i++) {
        FingerprintType f = fp[i + 2 * blockLength];
        if (f != 0) {
            bits[i >> 6] |= (1L << (i & 63));
            setBits++;
        }
    }
    fingerprints = new FingerprintType[2 * blockLength + setBits];
    for (size_t i = 0; i < 2 * blockLength; i++) {
        fingerprints[i] = fp[i];
    }
    // there are setBits uninitialized values?
    for (size_t i = 2 * blockLength; i < 2 * blockLength + setBits; i++) {
        fingerprints[i] = 0; // to be verified?
    }
    for (size_t i = 2 * blockLength, j = i; i < 3 * blockLength;) {
        FingerprintType f = fp[i++];
        if (f != 0) {
            fingerprints[j++] = f;
        }
    }
    delete [] fp;
    rank = new Rank9(bits, bitCount);
    delete [] bits;
    totalSizeInBytes = (2 * blockLength + setBits) * sizeof(FingerprintType)
        + rank->getBitCount() / 8;
    return Ok;
}

template <typename ItemType, typename FingerprintType,
          typename HashFamily>
Status XorFilterPlus<ItemType, FingerprintType, HashFamily>::Contain(
    const ItemType &key) const {
    uint64_t hash = (*hasher)(key);
    FingerprintType f = (FingerprintType) fingerprint(hash);
    uint32_t r0 = (uint32_t) hash;
    uint32_t r1 = (uint32_t) rotl64(hash, 21);
    uint32_t r2 = (uint32_t) rotl64(hash, 42);
    uint32_t h0 = reduce(r0, blockLength);
    uint32_t h1 = reduce(r1, blockLength) + blockLength;
    uint32_t h2a = reduce(r2, blockLength);
    f ^= fingerprints[h0] ^ fingerprints[h1];
    uint64_t AndPartialRank = rank->getAndPartialRank(h2a);
    if ((AndPartialRank & 1) == 1) {
        uint32_t h2x = (uint32_t) ((AndPartialRank >> 1) + rank->remainingRank(h2a));
        f ^= fingerprints[h2x + 2 * blockLength];
    }
    return f == 0 ? Ok : NotFound;
}

template <typename ItemType, typename FingerprintType,
          typename HashFamily>
std::string XorFilterPlus<ItemType, FingerprintType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "XorFilterPlus Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  return ss.str();
}
}  // namespace xorfilter_plus
#endif  // XOR_FILTER_PLUS_XOR_FILTER_PLUS_H_
