#ifndef VQ_FILTER_VQ_FILTER_H_
#define VQ_FILTER_VQ_FILTER_H_

#include <assert.h>
#include <algorithm>

#include "hashutil.h"

#include "vqf_filter.h"
#include "vqf_filter.c"

using namespace std;
using namespace hashing;

namespace vqfilter {
// status returned by a VQ filter operation
enum Status {
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

// #define SORTED_ADD

template <typename ItemType, typename HashFamily = SimpleMixSplit>
class VQFilter {

  vqf_filter *filter;
  uint64_t bytesUsed;
  uint64_t range;
  double bitsPerItem;
  HashFamily hasher;

  double BitsPerItem() const { return bitsPerItem; }
  
  void ApplyBlock(uint64_t *tmp, int block, int len);

 public:
  explicit VQFilter(const size_t n) : hasher() {

#ifdef SORTED_ADD
    // when inserting in sorted order
    uint64_t nslots = (uint64_t) (n / 0.85);
#else
    // when inserting in random order
    uint64_t nslots = (uint64_t) (n / 0.93);
#endif

    if ((filter = vqf_init(nslots)) == NULL) {
      std::cout << "Can't allocate.\n";
      abort();
    }
    range = filter->metadata.range;
    bytesUsed = filter->metadata.total_size_in_bytes;
    bitsPerItem = (double) bytesUsed / n;

  }

  ~VQFilter() {
      free(filter);
  }

  // Add an item to the filter.
  Status Add(const ItemType &item);
  
  Status AddAll(const vector<ItemType> &data, const size_t start, const size_t end) {
    return AddAll(data.data(), start, end);
  }

  // Add an item to the filter.
  Status AddAll(const ItemType *data, const size_t start, const size_t end);
  
  // Report if the item is inserted, with false positive rate.
  Status Contain(const ItemType &item) const;

  /* methods for providing stats  */
  // summary infomation
  std::string Info() const;

  // number of current inserted items;
  size_t Size() const { return 0; }

  // size of the filter in bytes.
  size_t SizeInBytes() const { return bytesUsed; }
};

template <typename ItemType, typename HashFamily>
Status VQFilter<ItemType, HashFamily>::Add(
    const ItemType &key) {
    uint64_t hash = hasher(key);
    bool ret = vqf_insert(filter, hash);
    if (!ret) {
        std::cout << "failed insertion for key.\n";
        abort();
    }
    return Ok;
}

template <typename ItemType, typename HashFamily>
Status VQFilter<ItemType, HashFamily>::Contain(
    const ItemType &key) const {
    uint64_t hash = hasher(key);
    bool ret = vqf_is_present(filter, hash);
    return ret ? Ok : NotFound;
}

const int blockShift = 15;
const int blockLen = 1 << blockShift;

template <typename ItemType, typename HashFamily>
void VQFilter<ItemType, HashFamily>::ApplyBlock(uint64_t *tmp, int block, int len) {
  // std::cout << "addAll ApplyBlock block " << block << " len " << len << "\n";
  for (int i = 0; i < len; i++) {
    uint64_t hash = tmp[(block << blockShift) + i];
    // std::cout << "inserting " << hash << "\n";
    bool ret = vqf_insert(filter, hash);
    if (!ret) {
        std::cout << "failed insertion for key.\n";
        abort();
    }
  }
}

template <typename ItemType, typename HashFamily>
Status VQFilter<ItemType, HashFamily>::AddAll(
    const ItemType* keys, const size_t start, const size_t end) {
#ifdef SORTED_ADD
  int blocks = 1 + (end - start) / blockLen;
  uint64_t *tmp = new uint64_t[blocks * blockLen];
  int *tmpLen = new int[blocks]();
  // std::cout << "addAll blocks " << blocks << "\n";
  for (size_t i = start; i < end; i++) {
    uint64_t key = keys[i];
    uint64_t hash = hasher(key);
    // __uint128_t x = (__uint128_t)key * (__uint128_t)blocks;
    __uint128_t x = (__uint128_t)hash * (__uint128_t)blocks;
    int block = (uint64_t)(x >> 64);
    int len = tmpLen[block];
    tmp[(block << blockShift) + len] = hash;
    tmpLen[block] = len + 1;
    if (len + 1 == blockLen) {
      ApplyBlock(tmp, block, len + 1);
      tmpLen[block] = 0;
    }
  }
  for (int block = 0; block < blocks; block++) {
    ApplyBlock(tmp, block, tmpLen[block]);
    tmpLen[block] = 0;
  }
  delete[] tmp;
  delete[] tmpLen;
#else
  // unsorted
  for (size_t i = start; i < end; i++) {
    uint64_t key = keys[i];
    uint64_t hash = hasher(key);
    bool ret = vqf_insert(filter, hash);
    if (!ret) {
      std::cout << "failed insertion for key.\n";
      abort();
    }
  }
#endif
  return Ok;
}

template <typename ItemType, typename HashFamily>
std::string VQFilter<ItemType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "VQFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  if (Size() > 0) {
    ss << "\t\tk:   " << BitsPerItem() << "\n";
  } else {
    ss << "\t\tk:   N/A\n";
  }
  return ss.str();
}
}  // namespace vqfilter
#endif  // VQ_FILTER_VQ_FILTER_H_
