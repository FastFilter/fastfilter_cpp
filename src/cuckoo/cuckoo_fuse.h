#ifndef CUCKOO_FILTER_CUCKOO_FUSE_FILTER_H_
#define CUCKOO_FILTER_CUCKOO_FUSE_FILTER_H_

#include <assert.h>
#include <algorithm>

#include "debug.h"
#include "hashutil.h"
#include "packedtable.h"
#include "printutil.h"
#include "singletable.h"

namespace cuckoofusefilter {
// status returned by a cuckoo filter operation
enum Status {
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

#define FINGERPRINT_TAG

// maximum number of cuckoo kicks before claiming failure
const size_t kMaxCuckooCount = 200000;

// 
const size_t segmentLengthBits = 14;
const size_t segmentLength = 1 << segmentLengthBits;

// A cuckoo filter class exposes a Bloomier filter interface,
// providing methods of Add, Delete, Contain.
template <typename ItemType, typename FingerprintType,
          typename HashFamily = hashing::SimpleMixSplit>
class CuckooFuseFilter {

  size_t size;
  size_t arrayLength;
  size_t segmentCount;
  FingerprintType *fingerprints;

  HashFamily* hasher;

  inline FingerprintType fingerprint(const uint64_t hash) const {
#ifdef FINGERPRINT_TAG
    return (FingerprintType) (hash >> 16);
#else
    FingerprintType x = (FingerprintType) (hash >> 16);
    x += (x == 0);    
    return x;
#endif
  }
  
  inline void segmentAndPos(uint64_t hash, int* seg, int* hh) const {
    __uint128_t x = (__uint128_t) hash * (__uint128_t) segmentCount;
    *seg = (uint64_t)(x >> 64);
    *hh = (hash >> 32) ^ hash;
  }

  inline uint64_t fingerprintHash(const FingerprintType fp) const {
    uint64_t x = fp;
#ifdef FINGERPRINT_TAG  
    x >>= 2;
#endif
    x *= 0x5bd1e995;
    return x ^ (x >> 32);
  }

 public:
 
  explicit CuckooFuseFilter(const size_t size) {
    hasher = new HashFamily();
    this->size = size;
    size_t capacity = size / 0.80 + 3 * segmentLength;
    this->segmentCount = (capacity - 3 * segmentLength) / segmentLength;
    this->arrayLength = (segmentCount + 3) * segmentLength;
    fingerprints = new FingerprintType[arrayLength]();
    std::fill_n(fingerprints, arrayLength, 0);
  }

  ~CuckooFuseFilter() {
    delete[] fingerprints;
    delete hasher;
  }

  // Add an item to the filter.
  Status Add(const ItemType &item);

  // Report if the item is inserted, with false positive rate.
  Status Contain(const ItemType &item) const;

  // Delete an key from the filter
  Status Delete(const ItemType &item);

  /* methods for providing stats  */
  // summary infomation
  std::string Info() const;

  // number of current inserted items;
  size_t Size() const { return size; }

  // size of the filter in bytes.
  size_t SizeInBytes() const { return arrayLength * sizeof(FingerprintType); }
};

template <typename ItemType, typename FingerprintType, typename HashFamily>
Status CuckooFuseFilter<ItemType, FingerprintType, HashFamily>::Add(
    const ItemType &key) {
    
    uint64_t hash = (*hasher)(key);
    int seg, hh;
    segmentAndPos(hash, &seg, &hh);
#ifdef FINGERPRINT_TAG
    FingerprintType fp = fingerprint(hash) << 2;
    uint64_t fh = fingerprintHash(fp);
    int h0 = (seg + 0) * segmentLength + (size_t)((hh) & (segmentLength - 1));
    int h1 = (seg + 1) * segmentLength + (size_t)((hh ^ fh) & (segmentLength - 1));
    int h2 = (seg + 2) * segmentLength + (size_t)((hh ^ (fh >> 32)) & (segmentLength - 1));
#else
    FingerprintType fp = fingerprint(hash);
    uint64_t fh = fingerprintHash(fp);
    hh = (hh << 2) ^ fh;
    int h0 = (seg + 0) * segmentLength + (size_t)((hh) & (segmentLength - 1));
    int h1 = (seg + 1) * segmentLength + (size_t)((hh ^ (fh << 2) ^ 1) & (segmentLength - 1));
    int h2 = (seg + 2) * segmentLength + (size_t)((hh ^ (fh >> 16 << 2) ^ 2) & (segmentLength - 1));
#endif
    for (uint32_t count = 0; count < kMaxCuckooCount; count++) {
#ifdef FINGERPRINT_TAG
        if (fingerprints[h0] == 0) {
            fingerprints[h0] = fp | 1;
            return Ok;
        } else if (fingerprints[h1] == 0) {
            fingerprints[h1] = fp | 2;
            return Ok;
        } else if (fingerprints[h2] == 0) {
            fingerprints[h2] = fp | 3;
            return Ok;
        }
#else
        if (fingerprints[h0] == 0) {
            fingerprints[h0] = fp;
            return Ok;
        } else if (fingerprints[h1] == 0) {
            fingerprints[h1] = fp;
            return Ok;
        } else if (fingerprints[h2] == 0) {
            fingerprints[h2] = fp;
            return Ok;
        }
#endif
        size_t m = rand() % 3;
        size_t old;
        if (m == 0) {
            old = h0;
        } else if (m == 1) {
            old = h1;
        } else {
            old = h2;
        }
#ifdef FINGERPRINT_TAG
        size_t m2 = fingerprints[old] & 0x3;
        seg = seg - (m2 - 1);
        FingerprintType fp2 = fingerprints[old] >> 2 << 2;
        if (m == 0) {
            fingerprints[h0] = fp | 1;
        } else if (m == 1) {
            fingerprints[h1] = fp | 2;
        } else if (m == 2) {
            fingerprints[h2] = fp | 3;
        }
        fp = fp2;
        uint64_t fh = fingerprintHash(fp);
        if (m2 == 1) {
            h0 = old;
            h1 = (h0 + (1 * segmentLength)) ^ (fh & (segmentLength - 1));
            h2 = (h0 + (2 * segmentLength)) ^ ((fh >> 32) & (segmentLength - 1));
        } else if (m2 == 2) {
            h1 = old;
            h0 = (h1 - (1 * segmentLength)) ^ (fh & (segmentLength - 1));
            h2 = (h0 + (2 * segmentLength)) ^ ((fh >> 32) & (segmentLength - 1));
        } else if (m2 == 3) {
            h2 = old;
            h0 = (h2 - (2 * segmentLength)) ^ ((fh >> 32) & (segmentLength - 1));
            h1 = (h0 + (1 * segmentLength)) ^ (fh & (segmentLength - 1));
        }
#else
        size_t m2 = (old ^ fingerprintHash(fingerprints[old])) & 0x3;
        
        seg = seg - m2;
        FingerprintType fp2 = fingerprints[old];
        if (m == 0) {
            fingerprints[h0] = fp;
        } else if (m == 1) {
            fingerprints[h1] = fp;
        } else if (m == 2) {
            fingerprints[h2] = fp;
        }
        fp = fp2;
        uint64_t fh = fingerprintHash(fp);
        if (m2 == 0) {
            h0 = old;
            h1 = (h0 + (1 * segmentLength)) ^ (((fh << 2) ^ 1) & (segmentLength - 1));
            h2 = (h0 + (2 * segmentLength)) ^ (((fh >> 16 << 2) ^ 2) & (segmentLength - 1));
        } else if (m2 == 1) {
            h1 = old;
            h0 = (h1 - (1 * segmentLength)) ^ (((fh << 2) ^ 1) & (segmentLength - 1));
            h2 = (h0 + (2 * segmentLength)) ^ (((fh >> 16 << 2) ^ 2) & (segmentLength - 1));
        } else if (m2 == 2) {
            h2 = old;
            h0 = (h2 - (2 * segmentLength)) ^ (((fh >> 16 << 2) ^ 2) & (segmentLength - 1));
            h1 = (h0 + (1 * segmentLength)) ^ (((fh << 2) ^ 1) & (segmentLength - 1));
        }
#endif
    }
    return NotEnoughSpace;
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
Status CuckooFuseFilter<ItemType, FingerprintType, HashFamily>::Contain(
    const ItemType &key) const {

    uint64_t hash = (*hasher)(key);
    int seg, hh;
    segmentAndPos(hash, &seg, &hh);
#ifdef FINGERPRINT_TAG
    FingerprintType fp = fingerprint(hash) << 2;
    uint64_t fh = fingerprintHash(fp);
    int h0 = (seg + 0) * segmentLength + (size_t)((hh) & (segmentLength - 1));
    int h1 = (seg + 1) * segmentLength + (size_t)((hh ^ fh) & (segmentLength - 1));
    int h2 = (seg + 2) * segmentLength + (size_t)((hh ^ (fh >> 32)) & (segmentLength - 1));
    FingerprintType f0 = fingerprints[h0];
    FingerprintType f1 = fingerprints[h1];
    FingerprintType f2 = fingerprints[h2];
    return ((f0 == (fp ^ 1)) | (f1 == (fp ^ 2)) | (f2 == (fp ^ 3))) ? Ok : NotFound;
#else
    FingerprintType fp = fingerprint(hash);
    uint64_t fh = fingerprintHash(fp);
    hh = (hh << 2) ^ fh;
    int h0 = (seg + 0) * segmentLength + (size_t)((hh) & (segmentLength - 1));
    int h1 = (seg + 1) * segmentLength + (size_t)((hh ^ (fh << 2) ^ 1) & (segmentLength - 1));
    int h2 = (seg + 2) * segmentLength + (size_t)((hh ^ (fh >> 16 << 2) ^ 2) & (segmentLength - 1));
    FingerprintType f0 = fingerprints[h0];
    FingerprintType f1 = fingerprints[h1];
    FingerprintType f2 = fingerprints[h2];
    return ((f0 == fp) | (f1 == fp) | (f2 == fp)) ? Ok : NotFound;
#endif
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
Status CuckooFuseFilter<ItemType, FingerprintType, HashFamily>::Delete(
    const ItemType &key) {
    
    uint64_t hash = (*hasher)(key);
    int seg, hh;
    segmentAndPos(hash, &seg, &hh);
#ifdef FINGERPRINT_TAG
    FingerprintType fp = fingerprint(hash) << 2;
    uint64_t fh = fingerprintHash(fp);
    int h0 = (seg + 0) * segmentLength + (size_t)((hh) & (segmentLength - 1));
    int h1 = (seg + 1) * segmentLength + (size_t)((hh ^ fh) & (segmentLength - 1));
    int h2 = (seg + 2) * segmentLength + (size_t)((hh ^ (fh >> 32)) & (segmentLength - 1));
    if (fingerprints[h0] == (fp | 1)) {
        fingerprints[h0] = 0;
        return Ok;
    } else if (fingerprints[h1] == (fp | 2)) {
        fingerprints[h1] = 0;
        return Ok;
    } else if (fingerprints[h2] == (fp | 3)) {
        fingerprints[h2] = 0;
        return Ok;
    }
#else
    FingerprintType fp = fingerprint(hash);
    uint64_t fh = fingerprintHash(fp);
    hh = (hh << 2) ^ fh;
    int h0 = (seg + 0) * segmentLength + (size_t)((hh) & (segmentLength - 1));
    int h1 = (seg + 1) * segmentLength + (size_t)((hh ^ (fh << 2) ^ 1) & (segmentLength - 1));
    int h2 = (seg + 2) * segmentLength + (size_t)((hh ^ (fh >> 16 << 2) ^ 2) & (segmentLength - 1));
    if (fingerprints[h0] == fp) {
        fingerprints[h0] = 0;
        return Ok;
    } else if (fingerprints[h1] == fp) {
        fingerprints[h1] = 0;
        return Ok;
    } else if (fingerprints[h2] == fp) {
        fingerprints[h2] = 0;
        return Ok;
    }
#endif
    return NotFound;
}

template <typename ItemType, typename FingerprintType, typename HashFamily>
std::string CuckooFuseFilter<ItemType, FingerprintType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "CuckooFuseFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n"
     << "\t\tHashtable size: " << SizeInBytes() << " KB\n";
  return ss.str();
}
}  // namespace cuckoofusefilter
#endif  // CUCKOO_FILTER_CUCKOO_FUSE_FILTER_H_
