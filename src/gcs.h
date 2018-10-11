#ifndef GCS_FILTER_GCS_FILTER_H_
#define GCS_FILTER_GCS_FILTER_H_

#include <assert.h>
#include <algorithm>

#include "debug.h"
#include "hashutil.h"
#include "printutil.h"

using namespace std;
using namespace cuckoofilter;

namespace gcsfilter {
// status returned by a gcs filter operation
enum Status {
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

inline uint32_t fingerprint(uint64_t hash) {
    return (uint32_t) (hash & ((1 << 8) - 1));
}

inline uint32_t reduce(uint32_t hash, uint32_t n) {
    // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
    return (uint32_t) (((uint64_t) hash * n) >> 32);
}

// MultiStageMonotoneList

inline int numberOfLeadingZeros64(uint64_t x) {
    // If x is 0, the result is undefined.
    return x == 0 ? 64 : __builtin_clzl(x);
}

inline int numberOfLeadingZeros32(uint32_t x) {
    // If x is 0, the result is undefined.
    return x == 0 ? 32 : __builtin_clz(x);
}

inline uint64_t readNumber(uint64_t* data, uint64_t pos, int bitCount) {
    if (bitCount == 0) {
        return 0;
    }
    int remainingBits = 64 - (pos & 63);
    int index = pos >> 6;
    long x = data[index];
    if (bitCount <= remainingBits) {
        x >>= remainingBits - bitCount;
        return x & ((1L << bitCount) - 1);
    }
    x = x & ((1L << remainingBits) - 1);
    return (x << (bitCount - remainingBits)) |
        (data[index + 1] >> (64 - bitCount + remainingBits));
}

uint64_t readUntilZeroMore(uint64_t* data, int count, int index) {
    while (true) {
        uint64_t x = data[++index];
        if (x == UINT64_MAX) {
            count += 64;
            continue;
        }
        return count + numberOfLeadingZeros64(~x);
    }
}

inline uint64_t readUntilZero(uint64_t* data, uint64_t pos) {
    int remainingBits = 64 - (pos & 63);
    int index = pos >> 6;
    uint64_t x = data[index] << (64 - remainingBits);
    int count = numberOfLeadingZeros64(~x);
    if (count < remainingBits) {
        return count;
    }
    return readUntilZeroMore(data, count, index);
}

uint64_t writeNumber(uint64_t* data, uint64_t pos, uint64_t x, int bitCount) {
    // while (bitCount-- > 0) {
    //     writeBit((x >>> bitCount) & 1);
    // }
    if (bitCount == 0) {
        return pos;
    }
    int remainingBits = 64 - (pos & 63);
    int index = pos >> 6;
    if (bitCount <= remainingBits) {
        data[index] |= x << (remainingBits - bitCount);
    } else {
        data[index] |= x >> (bitCount - remainingBits);
        data[index + 1] |= x << (64 - bitCount + remainingBits);
    }
    pos += bitCount;
    return pos;
}

uint64_t writeGolombRice(uint64_t* data, uint64_t pos, int shift, uint64_t value) {
    uint64_t q = value >> shift;
    assert(q < 63);
    uint64_t m = (2L << q) - 2;
    pos = writeNumber(data, pos, m, (int) (q + 1));
    pos = writeNumber(data, pos, value & ((1L << shift) - 1), shift);
    return pos;
}

uint64_t getScaleFactor(int64_t multiply, int64_t divide) {
    return divide == 0 ? 0 : ((uint64_t) multiply << 32) / divide + 1;
}

typedef struct MultiStageMonotoneList {
    uint64_t* data;
    uint32_t dataBits;
    uint64_t startLevel1, startLevel2, startLevel3;
    int bitCount1, bitCount2, bitCount3;
    uint32_t count1, count2, count3;
    uint64_t factor;
    int32_t add;
} MultiStageMonotoneList;

MultiStageMonotoneList list;

#define SHIFT1 6
#define SHIFT2 3
#define FACTOR1 32
#define FACTOR2 16
#define MAX_VALUE32 0x7fffffff

void MultiStageMonotoneList_generate(MultiStageMonotoneList* list, uint32_t* data, int length) {
    int count3 = length;
    // verify it is monotone
    for (int i = 1; i < count3; i++) {
        assert (data[i - 1] <= data[i]);
    }
    uint64_t diff = data[count3 - 1] - data[0];
    uint64_t factor = getScaleFactor(diff, count3);
    list->factor = factor;
    int32_t add = data[0];
    for (int i = 1; i < count3; i++) {
        int expected = (int) ((i * factor) >> 32);
        int x = data[i];
        add = min(add, x - expected);
    }
    list->count3 = count3;
    list->add = add;
    int count2 = (count3 + (1 << SHIFT2) - 1) >> SHIFT2;
    int count1 = (count3 + (1 << SHIFT1) - 1) >> SHIFT1;
    int* group1 = new int[count1];
    int* group2 = new int[count2];
    int* group3 = new int[count3];
    for (int i = 0; i < count3; i++) {
        // int expected = (int) (i * max / count3);
        int expected = (int) ((i * factor) >> 32) + add;
        int got = data[i];
        int x = got - expected;
        assert(x >= 0);
        group3[i] = x;
    }
    int a = MAX_VALUE32;
    for (int i = 0; i < count3; i++) {
        int x = group3[i];
        a = min(a, x);
        if ((i +1) >> SHIFT2 != i >> SHIFT2 || i == count3 - 1) {
            group2[i >> SHIFT2] = a / FACTOR2;
            a = MAX_VALUE32;
        }
    }
    a = MAX_VALUE32;
    for (int i = 0; i < count3; i++) {
        int d = group2[i >> SHIFT2] * FACTOR2;
        int x = group3[i];
        group3[i] -= d;
        assert(group3[i] >= 0);
        a = min(a, x);
        if ((i + 1) >> SHIFT1 != i >> SHIFT1 || i == count3 - 1) {
            group1[i >> SHIFT1] = a / FACTOR1;
            a = MAX_VALUE32;
        }
    }
    int last = -1;
    for (int i = 0; i < count3; i++) {
        int i2 = i >> SHIFT2;
        if (i2 == last) {
            continue;
        }
        int d = group1[i >> SHIFT1] * FACTOR1;
        group2[i2] -= d / FACTOR2;
        last = i2;
    }
    int max1 = 0, max2 = 0, max3 = 0;
    for (int i = 0; i < count3; i++) {
        max3 = max(max3, group3[i]);
    }
    for (int i = 0; i < count2; i++) {
        max2 = max(max2, group2[i]);
    }
    for (int i = 0; i < count1; i++) {
        max1 = max(max1, group1[i]);
    }
    int bitCount1 = 32 - numberOfLeadingZeros32(max1);
    int bitCount2 = 32 - numberOfLeadingZeros32(max2);
    int bitCount3 = 32 - numberOfLeadingZeros32(max3);
    list->bitCount1 = bitCount1;
    list->bitCount2 = bitCount2;
    list->bitCount3 = bitCount3;
    int pos = 0;
    list->dataBits = bitCount1 * count1 + bitCount2 * count2 + bitCount3 * count3;
    list->startLevel1 = pos;
    size_t wordlen = (list->dataBits + 63) / 64;
    list->data = new uint64_t[wordlen];
    memset(list->data,0,wordlen*sizeof(uint64_t)); // presumably this is needed
    for (int i = 0; i < count1; i++) {
        pos = writeNumber(list->data, pos, group1[i], bitCount1);
    }
    list->startLevel2 = pos;
    for (int i = 0; i < count2; i++) {
        pos = writeNumber(list->data, pos, group2[i], bitCount2);
    }
    list->startLevel3 = pos;
    for (int i = 0; i < count3; i++) {
        pos = writeNumber(list->data, pos, group3[i], bitCount3);
    }
    delete[] group1;
    delete[] group2;
    delete[] group3;
}

inline uint32_t MultiStageMonotoneList_get(const MultiStageMonotoneList* list, uint32_t i) {
    int expected = (int) ((i * list->factor) >> 32) + list->add;
    long a = readNumber(list->data, list->startLevel1 + (i >> SHIFT1) * list->bitCount1, list->bitCount1);
    long b = readNumber(list->data, list->startLevel2 + (i >> SHIFT2) * list->bitCount2, list->bitCount2);
    long c = readNumber(list->data, list->startLevel3 + i * list->bitCount3, list->bitCount3);
    return (int) (expected + a * FACTOR1 + b * FACTOR2 + c);
}

template <typename ItemType, size_t bits_per_item,
          typename HashFamily = TwoIndependentMultiplyShift>
class GcsFilter {

  int golombShift;
  int bufferSize;
  int bucketCount;
  int fingerprintMask;
  MultiStageMonotoneList monotoneList;
  int startBuckets;
  uint64_t* bucketData;
  size_t bucketDataBits;

  HashFamily hasher;

  double BitsPerItem() const { return 8.0; }

 public:
  explicit GcsFilter(const size_t len) : hasher() {
  }

  ~GcsFilter() {
    delete[] bucketData;
    delete[] monotoneList.data;
  }

  Status AddAll(const vector<ItemType> data, const size_t start, const size_t end);

  // Report if the item is inserted, with false positive rate.
  Status Contain(const ItemType &item) const;

  /* methods for providing stats  */
  // summary infomation
  std::string Info() const;

  // number of current inserted items;
  size_t Size() const { return 0; }

  // size of the filter in bytes.
  size_t SizeInBytes() const { return monotoneList.dataBits / 8 + bucketDataBits / 8; }
};

int compare_uint64(const void* a, const void* b) {
    const uint64_t ai = *(const uint64_t*)a;
    const uint64_t bi = *(const uint64_t*)b;
    return ai < bi ? -1 : ai > bi ? 1 : 0;
}

template <typename ItemType, size_t bits_per_item,
          typename HashFamily>
Status GcsFilter<ItemType, bits_per_item, HashFamily>::AddAll(
    const vector<ItemType> keys, const size_t start, const size_t end) {

    int len = end - start;
    // this was found experimentally
    int fingerprintBits = bits_per_item;
    golombShift = fingerprintBits - 1;
    int averageBucketSize = 16;
    // due to average bucket size of 16
    fingerprintBits += 4;
    uint64_t* data = new uint64_t[len];
    fingerprintMask = (1 << fingerprintBits) - 1;
    bucketCount = (int) ((len + averageBucketSize - 1) / averageBucketSize);
    for (int i = 0; i < len; i++) {
        uint64_t h = hasher(keys[i + start]);
        uint64_t b = reduce((int) (h >> 32), bucketCount);
        data[i] = (b << 32) | (h & fingerprintMask);
    }
    qsort(data, len, sizeof(uint64_t), compare_uint64);
    size_t bucketslen = 10L * fingerprintBits * len / 64;
    uint64_t* buckets = new uint64_t[bucketslen];
    memset(buckets, 0, sizeof(uint64_t[bucketslen]));
    uint32_t* startList = new uint32_t[bucketCount + 1];
    memset(startList, 0, sizeof(uint32_t[bucketCount + 1]));
    int bucket = 0;
    long last = 0;
    int pos = 0;
    for (int i = 0; i < len; i++) {
        long x = data[i];
        int b = (int) (x >> 32);
        while (bucket <= b) {
            startList[bucket++] = pos;
            last = 0;
        }
        x &= fingerprintMask;
        long diff = x - last;
        last = x;
        pos = writeGolombRice(buckets, pos, golombShift, diff);
    }
    while (bucket <= bucketCount) {
        startList[bucket++] = pos;
    }

    this->bucketData = buckets;
    startBuckets = 0;
    bucketDataBits = pos;
    MultiStageMonotoneList_generate(&monotoneList, startList, bucketCount + 1);

    delete[] data;
    delete[] startList;

    return Ok;
}

template <typename ItemType, size_t bits_per_item,
          typename HashFamily>
Status GcsFilter<ItemType, bits_per_item, HashFamily>::Contain(
    const ItemType &key) const {
    uint64_t hashCode = hasher(key);
    size_t b = reduce((uint32_t) (hashCode >> 32), bucketCount);
    uint64_t fingerprint = hashCode & fingerprintMask;

    size_t p = MultiStageMonotoneList_get(&monotoneList, b);
    size_t startNext = MultiStageMonotoneList_get(&monotoneList, b + 1);

    uint64_t x = 0;
    while (p < startNext) {
        long q = readUntilZero(bucketData, p);
        p += q + 1;
        x += (q << golombShift) | readNumber(bucketData, p, golombShift);
        if (x == fingerprint) {
            return Ok;
        } else if (x > fingerprint) {
            break;
        }
        p += golombShift;
    }
    return NotFound;
}

template <typename ItemType, size_t bits_per_item,
          typename HashFamily>
std::string GcsFilter<ItemType, bits_per_item, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "GcsFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  if (Size() > 0) {
    ss << "\t\tbit/key:   " << BitsPerItem() << "\n";
  } else {
    ss << "\t\tbit/key:   N/A\n";
  }
  return ss.str();
}
}  // namespace gcsfilter
#endif  // GCS_FILTER_GCS_FILTER_H_
