#ifndef XOR_FILTER_PLUS_XOR_FILTER_PLUS_H_
#define XOR_FILTER_PLUS_XOR_FILTER_PLUS_H_

#include <assert.h>
#include <algorithm>

#include "hashutil.h"

using namespace std;
using namespace hashing;

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

// from https://github.com/rob-p/rank_speed_test/blob/master/src/poppy.cpp
// license:
// GNU Lesser General Public License v3.0
// which looks like it's from
// https://github.com/efficient/rankselect
// license:
// Copyright (C) 2013, Carnegie Mellon University
// Licensed under the Apache License, Version 2.0 (the "License")

const int kWordSize = 64;
const int kBasicBlockSize = 512;
const int kBasicBlockBits = 9;
const int kBasicBlockMask = kBasicBlockSize - 1;
const int kWordCountPerBasicBlock = kBasicBlockSize / kWordSize;
const int kCacheLineSize = 64;

// #define USE_POPPY

class Poppy {
private:
    uint64_t *bits_;
    uint64_t num_bits_;
    uint64_t num_counts_;

    uint64_t *l2Entries_;
    uint64_t l2EntryCount_;
    uint64_t *l1Entries_;
    uint64_t l1EntryCount_;
    uint64_t basicBlockCount_;

    uint32_t *loc_[1 << 16];
    uint32_t locCount_[1 << 16];

    static const int kLocFreq = 8192;
    static const int kLocFreqMask = 8191;
    static const int kL2EntryCountPerL1Entry = 1 << 21;
    static const int kBasicBlockCountPerL1Entry = 1 << 23;

public:
    Poppy(uint64_t * const bits, const uint64_t num_bits);

    inline uint64_t rank(uint64_t pos);
    inline uint64_t get(uint64_t pos) {
        return (bits_[(size_t) (pos >> 6)] >> (63 - (pos & 63))) & 1;
    }

    uint64_t getBitCount() {
        return num_counts_;
    }

};

#define popcountsize 64ULL
#define popcountmask (popcountsize - 1)
#define _mm_popcnt_u64 bitCount64

inline uint64_t popcountLinear(uint64_t *bits, uint64_t x, uint64_t nbits) {
    if (nbits == 0)
        return 0;

    uint64_t lastword = (nbits - 1) / popcountsize;
    uint64_t p = 0;

    for (int i = 0; i < lastword; i++) {
        p += _mm_popcnt_u64(bits[x+i]);
    }

    uint64_t lastshifted = bits[x+lastword] >> (63 - ((nbits - 1) & popcountmask));
    p += _mm_popcnt_u64(lastshifted);
    return p;
}

Poppy::Poppy(uint64_t * const bits, uint64_t num_bits) {
    size_t bitsArraySize = (size_t) ((num_bits + 511) / 512) * 512 / 64;
    // TODO use posix_memalign
    posix_memalign((void **) &bits_, kCacheLineSize, (bitsArraySize + 16) * sizeof(uint64_t));
    // bits_ = new uint64_t[bitsArraySize + 16]();

    for(int i=0; i<bitsArraySize; i++) {
        uint64_t x = bits[i];
        uint64_t y = 0;
        for(int j=0; j<64; j++) {
            y = (y << 1) | (x & 0x1);
            x >>= 1;
        }
        bits_[i] = y;
    }
    num_bits = (bitsArraySize + 16) * 64;

    // memcpy(bits_, bits, (bitsArraySize - 1) * sizeof(uint64_t));
    // bits_[bitsArraySize - 1] = 0;

//    bits_ = bits;
// std::cout << "poppy.init0 " << num_bits << "\n";


    num_bits_ = num_bits;
    num_counts_ = 0;

    l1EntryCount_ = std::max(num_bits_ >> 32, (uint64_t) 1);
    l2EntryCount_ = num_bits_ >> 11;
    basicBlockCount_ = num_bits_ / kBasicBlockSize;

    // assert(
    posix_memalign((void **) &l1Entries_, kCacheLineSize, l1EntryCount_ * sizeof(uint64_t));
    //  >= 0);
    // assert(
    posix_memalign((void **) &l2Entries_, kCacheLineSize, l2EntryCount_ * sizeof(uint64_t));
    //  >= 0);

    uint64_t l2Id = 0;
    uint64_t basicBlockId = 0;

    memset(locCount_, 0, sizeof(locCount_));

    for (uint64_t i = 0; i < l1EntryCount_; i++) {
        l1Entries_[i] = num_counts_;
        uint32_t cum = 0;
        for (int k = 0; k < kL2EntryCountPerL1Entry; k++) {
            l2Entries_[l2Id] = cum;
            for (int offset = 0; offset < 30; offset += 10) {
                int c = popcountLinear(bits_,
                                       basicBlockId * kWordCountPerBasicBlock,
                                       kBasicBlockSize);
                cum += c;
                basicBlockId++;
                l2Entries_[l2Id] |= (uint64_t) c << (32 + offset);
            }
            cum += popcountLinear(bits_, basicBlockId * kWordCountPerBasicBlock, kBasicBlockSize);
            basicBlockId++;

            if (++l2Id >= l2EntryCount_) break;
        }

        locCount_[i] = (cum + kLocFreq - 1) / kLocFreq;
        num_counts_ += cum;
    }
    basicBlockId = 0;
    for (uint64_t i = 0; i < l1EntryCount_; i++) {
        loc_[i] = new uint32_t[locCount_[i]];
        locCount_[i] = 0;

        uint32_t oneCount = 0;

        for (uint32_t k = 0; k < kBasicBlockCountPerL1Entry; k++) {
            uint64_t woff = basicBlockId * kWordCountPerBasicBlock;
            for (int widx = 0; widx < kWordCountPerBasicBlock; widx++)
                for (int bit = 0; bit < kWordSize; bit++)
                    if (bits_[woff + widx] & (1ULL << (63 - bit))) {
                        oneCount++;
                        if ((oneCount & kLocFreqMask) == 1) {
                            loc_[i][locCount_[i]] = k * kBasicBlockSize + widx * kWordSize + bit;
                            locCount_[i]++;
                        }
                    }

            basicBlockId++;
            if (basicBlockId >= basicBlockCount_) break;
        }
    }

// to ensure everything is OK
Rank9 *r = new Rank9(bits, num_bits);
for(int i=0; i<num_bits; i++) {
    if (r->rank(i) != rank(i)) {
        std::cout << "rank " << i << " of " << num_bits << " r9 " << r->rank(i) << " poppy " << rank(i) << "\n";
        break;
    }
}
for(int i=0; i<num_bits; i++) {
    if (r->get(i) != get(i)) {
        std::cout << "get " << i << " of " << num_bits << " r9 " << r->get(i) << " poppy " << get(i) << "\n";
        break;
    }
}
delete r;

}

inline uint64_t Poppy::rank(uint64_t pos) {
    // assert(pos <= num_bits_);
//    --pos;
//        std::cout << "poppy.rank " << pos << "\n";


    uint64_t l1Id = pos >> 32;
    uint64_t l2Id = pos >> 11;
    uint64_t x = l2Entries_[l2Id];

    uint64_t res = l1Entries_[l1Id] + (x & 0xFFFFFFFFULL);
    x >>= 32;

    int groupId = (pos & 2047) / 512;
    for (int i = 0; i < groupId; i++) {
        res += x & 1023;
        x >>= 10;
    }
    res += popcountLinear(bits_, (l2Id * 4 + groupId) * kWordCountPerBasicBlock, (pos & 511));

    return res;
}

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
#ifdef USE_POPPY
  Poppy *rank = NULL;
#else
  Rank9 *rank = NULL;
#endif
  size_t totalSizeInBytes;

  HashFamily* hasher;

  inline FingerprintType fingerprint(const uint64_t hash) const {
    return (FingerprintType) (hash ^ (hash >> 32));
  }

 public:
  explicit XorFilterPlus(const size_t size) {
    hasher = new HashFamily();
    this->size = size;
    this->arrayLength = 32 + 1.23 * size;
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

  Status AddAll(const vector<ItemType>& data, const size_t start, const size_t end) {
      return AddAll(data.data(), start, end);
  }
  Status AddAll(const ItemType * data, const size_t start, const size_t end);

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
    const ItemType* keys, const size_t start, const size_t end) {
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
            std::cout << (end - start) << " keys; arrayLength " << arrayLength
                << " blockLength " << blockLength
                << " reverseOrderPos " << reverseOrderPos << "\n";
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
    uint64_t *bits = new uint64_t[(bitCount + 63) / 64]();
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
#ifdef USE_POPPY
    rank = new Poppy(bits, bitCount);
#else
    rank = new Rank9(bits, bitCount);
#endif
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
#ifdef USE_POPPY
    if (rank->get(h2a)) {
        uint32_t h2x = (uint32_t) rank->rank(h2a);
        f ^= fingerprints[h2x + 2 * blockLength];
    }
#else
    uint64_t bitAndPartialRank = rank->getAndPartialRank(h2a);
    if ((bitAndPartialRank & 1) == 1) {
        uint32_t h2x = (uint32_t) ((bitAndPartialRank >> 1) + rank->remainingRank(h2a));
        f ^= fingerprints[h2x + 2 * blockLength];
    }
#endif
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
