#ifndef GQ_FILTER_GQ_FILTER_H_
#define GQ_FILTER_GQ_FILTER_H_

#include <assert.h>
#include <algorithm>

#include "debug.h"
#include "hashutil.h"
#include "printutil.h"

#include "gqf_hashutil.h"
#include "gqf_hashutil.c"
#include "gqf.h"
#include "gqf_int.h"
#include "gqf.c"

using namespace std;
using namespace cuckoofilter;

namespace gqfilter {
// status returned by a GQ filter operation
enum Status {
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

template <typename ItemType, size_t bits_per_item,
    typename HashFamily = TwoIndependentMultiplyShift>
class GQFilter {

  QF qf;
  uint64_t mask;
  uint64_t bytesUsed;
  double bitsPerItem;
  HashFamily hasher;

  double BitsPerItem() const { return bitsPerItem; }

 public:
  explicit GQFilter(const size_t n) : hasher() {

    uint64_t qbits;
    uint64_t nslots;
    qbits = 0;
    do {
        qbits++;
        nslots = (1ULL << qbits);
    } while(nslots * 0.9 < n);
    uint64_t nhashbits = qbits + 8;
    mask = (1ULL << nhashbits) - 1;

// std::cout << "(CQF: nslots " << nslots << " nhashbits " << nhashbits << " n " << n << " bitsPerItem " << bitsPerItem << ")\n";

// if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {
//    if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_DEFAULT, 0)) {
    if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_NONE, 0)) {
        std::cout << "Can't allocate CQF.\n";
        abort();
    }

    // qf_set_auto_resize(&qf, true);

    bytesUsed = qf.metadata->total_size_in_bytes;
    bitsPerItem = (double) bytesUsed / n;

  }

  ~GQFilter() {
      qf_free(&qf);
  }

  // Add an item to the filter.
  Status Add(const ItemType &item);

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

template <typename ItemType, size_t bits_per_item,
    typename HashFamily>
Status GQFilter<ItemType, bits_per_item, HashFamily>::Add(
    const ItemType &key) {
    uint64_t hash = hasher(key);
    // uint64_t hash = key;
    // int ret = qf_insert(&qf, hash & mask, 0, 1, QF_NO_LOCK | QF_KEY_IS_HASH);
    int ret = qf_insert(&qf, hash & mask, 0, 1, QF_NO_LOCK);
    if (ret < 0) {
        std::cout << "failed insertion for key.\n";
        if (ret == QF_NO_SPACE) {
            std::cout << "CQF is full.\n";
        } else if (ret == QF_COULDNT_LOCK) {
            std::cout << "TRY_ONCE_LOCK failed.\n";
        } else {
            std::cout << "Does not recognise return value.\n";
        }
        abort();
    }
    return Ok;
}

template <typename ItemType, size_t bits_per_item,
    typename HashFamily>
Status GQFilter<ItemType, bits_per_item, HashFamily>::Contain(
    const ItemType &key) const {
    uint64_t hash = hasher(key);
    // uint64_t hash = key;
    // uint64_t count = qf_count_key_value(&qf, hash & mask, 0, QF_NO_LOCK | QF_KEY_IS_HASH);
    uint64_t count = qf_count_key_value(&qf, hash & mask, 0, 0);
    return count > 0 ? Ok : NotFound;
}

template <typename ItemType, size_t bits_per_item,
    typename HashFamily>
std::string GQFilter<ItemType, bits_per_item, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "GQFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  if (Size() > 0) {
    ss << "\t\tk:   " << BitsPerItem() << "\n";
  } else {
    ss << "\t\tk:   N/A\n";
  }
  return ss.str();
}
}  // namespace gqfilter
#endif  // GQ_FILTER_GQ_FILTER_H_
