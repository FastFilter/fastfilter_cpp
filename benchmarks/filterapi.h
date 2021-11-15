#ifndef FILTERAPI_H
#define FILTERAPI_H
#include <climits>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <vector>
#include <set>
#include <stdio.h>

// morton
#include "compressed_cuckoo_filter.h"
#include "morton_sample_configs.h"
#include "cuckoofilter.h"
#include "cuckoofilter_stable.h"
#include "cuckoo_fuse.h"
#include "xorfilter.h"
#include "xorfilter_plus.h"
#include "xorfilter_singleheader.h"
#include "binaryfusefilter_singleheader.h"
#include "xor_binary_fuse_filter.h"
#include "bloom.h"
#include "counting_bloom.h"
#include "gcs.h"
#ifdef __AVX2__
#include "gqf_cpp.h"
#include "vqf_cpp.h"
#include "simd-block.h"
#endif
#include "simd-block-fixed-fpp.h"
#include "ribbon_impl.h"

using namespace std;
using namespace hashing;
using namespace cuckoofilter;
using namespace cuckoofusefilter;
using namespace xorfilter;
using namespace xorfilter_plus;
using namespace bloomfilter;
using namespace counting_bloomfilter;
using namespace gcsfilter;
using namespace CompressedCuckoo; // Morton filter namespace
#ifdef __AVX2__
using namespace gqfilter;
using namespace vqfilter;
#endif
using namespace ribbon;



// Inlining the "contains" which are executed within a tight loop can be both
// very detrimental or very beneficial, and which ways it goes depends on the
// compiler. It is unclear whether we want to benchmark the inlining of Contains,
// as it depends very much on how "contains" is used. So it is maybe reasonable
// to benchmark it without inlining.
//
#define CONTAIN_ATTRIBUTES  __attribute__ ((noinline))

template<typename Table>
struct FilterAPI {};

template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
struct FilterAPI<CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>> {
  using Table = CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table * table) {
    if (0 != table->Add(key)) {
      throw logic_error("The filter is too small to hold all of the elements");
    }
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Delete(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
struct FilterAPI<CuckooFilterStable<ItemType, bits_per_item, TableType, HashFamily>> {
  using Table = CuckooFilterStable<ItemType, bits_per_item, TableType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table * table) {
    if (0 != table->Add(key)) {
      throw logic_error("The filter is too small to hold all of the elements");
    }
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Delete(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<CuckooFuseFilter<ItemType, FingerprintType>> {
  using Table = CuckooFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table * table) {
    if (0 != table->Add(key)) {
      throw logic_error("The filter is too small to hold all of the elements");
    }
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Delete(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

#ifdef __aarch64__
template <typename HashFamily>
struct FilterAPI<SimdBlockFilterFixed<HashFamily>> {
  using Table = SimdBlockFilterFixed<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(add_count);
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};

#endif

#ifdef __AVX2__
template <typename HashFamily>
struct FilterAPI<SimdBlockFilter<HashFamily>> {
  using Table = SimdBlockFilter<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(ceil(log2(add_count)));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};

template <typename HashFamily>
struct FilterAPI<SimdBlockFilterFixed64<HashFamily>> {
  using Table = SimdBlockFilterFixed64<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(add_count);
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};

template <typename HashFamily>
struct FilterAPI<SimdBlockFilterFixed<HashFamily>> {
  using Table = SimdBlockFilterFixed<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(add_count);
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};
#endif

#ifdef __SSE41__
template <typename HashFamily>
struct FilterAPI<SimdBlockFilterFixed16<HashFamily>> {
  using Table = SimdBlockFilterFixed16<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(add_count);
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};
#endif

template <typename ItemType, typename FingerprintType>
struct FilterAPI<XorFilter<ItemType, FingerprintType>> {
  using Table = XorFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename CoeffType, bool kHomog, uint32_t kNumColumns, bool kSmash = false>
struct RibbonTS {
  static constexpr bool kIsFilter = true;
  static constexpr bool kHomogeneous = kHomog;
  static constexpr bool kFirstCoeffAlwaysOne = true;
  static constexpr bool kUseSmash = kSmash;
  using CoeffRow = CoeffType;
  using Hash = uint64_t;
  using Key = uint64_t;
  using Seed = uint32_t;
  using Index = size_t;
  using ResultRow = uint32_t;
  static constexpr bool kAllowZeroStarts = false;
  static constexpr uint32_t kFixedNumColumns = kNumColumns;

  static Hash HashFn(const Hash& input, Seed raw_seed) {
    // return input;
    uint64_t h = input + raw_seed;
    h ^= h >> 33;
    h *= UINT64_C(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h *= UINT64_C(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    return h;
  }
};

template <typename CoeffType, uint32_t kNumColumns, uint32_t kMilliBitsPerKey = 7700>
class HomogRibbonFilter {
  using TS = RibbonTS<CoeffType, /*kHomog*/ true, kNumColumns>;
  IMPORT_RIBBON_IMPL_TYPES(TS);

  size_t num_slots;
  size_t bytes;
  unique_ptr<char[]> ptr;
  InterleavedSoln soln;
  Hasher hasher;
public:
  static constexpr double kFractionalCols =
    kNumColumns == 0 ? kMilliBitsPerKey / 1000.0 : kNumColumns;

  static double GetBestOverheadFactor() {
    double overhead = (4.0 + kFractionalCols * 0.25) / (8.0 * sizeof(CoeffType));
    return 1.0 + overhead;
  }

  HomogRibbonFilter(size_t add_count)
      : num_slots(InterleavedSoln::RoundUpNumSlots((size_t)(GetBestOverheadFactor() * add_count))),
        bytes(static_cast<size_t>((num_slots * kFractionalCols + 7) / 8)),
        ptr(new char[bytes]),
        soln(ptr.get(), bytes) {
        }

  void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end) {
    Banding b(num_slots);
    (void)b.AddRange(keys.begin() + start, keys.begin() + end);
    soln.BackSubstFrom(b);
  }
  bool Contain(uint64_t key) const {
    return soln.FilterQuery(key, hasher);
  }
  size_t SizeInBytes() const {
    return bytes;
  }
};

template <typename CoeffType, uint32_t kNumColumns, uint32_t kMilliBitsPerKey>
struct FilterAPI<HomogRibbonFilter<CoeffType, kNumColumns, kMilliBitsPerKey>> {
  using Table = HomogRibbonFilter<CoeffType, kNumColumns, kMilliBitsPerKey>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Contain(key);
  }
};

template <typename CoeffType, uint32_t kNumColumns, uint32_t kMinPctOverhead, uint32_t kMilliBitsPerKey = 7700>
class BalancedRibbonFilter {
  using TS = RibbonTS<CoeffType, /*kHomog*/ false, kNumColumns>;
  IMPORT_RIBBON_IMPL_TYPES(TS);
  static constexpr uint32_t kBitsPerVshard = 8;
  using BalancedBanding = ribbon::BalancedBanding<TS, kBitsPerVshard>;
  using BalancedHasher = ribbon::BalancedHasher<TS, kBitsPerVshard>;

  uint32_t log2_vshards;
  size_t num_slots;

  size_t bytes;
  unique_ptr<char[]> ptr;
  InterleavedSoln soln;

  size_t meta_bytes;
  unique_ptr<char[]> meta_ptr;
  BalancedHasher hasher;
public:
  static constexpr double kFractionalCols =
    kNumColumns == 0 ? kMilliBitsPerKey / 1000.0 : kNumColumns;

  static double GetNumSlots(size_t add_count, uint32_t log2_vshards) {
    size_t add_per_vshard = add_count >> log2_vshards;

    double overhead;
    if (sizeof(CoeffType) == 8) {
      overhead = 0.0000055 * add_per_vshard; // FIXME?
    } else if (sizeof(CoeffType) == 4) {
      overhead = 0.00005 * add_per_vshard;
    } else if (sizeof(CoeffType) == 2) {
      overhead = 0.00010 * add_per_vshard; // FIXME?
    } else {
      assert(sizeof(CoeffType) == 16);
      overhead = 0.0000013 * add_per_vshard;
    }
    overhead = std::max(overhead, 0.01 * kMinPctOverhead);
    return InterleavedSoln::RoundUpNumSlots((size_t)(add_count + overhead * add_count + add_per_vshard / 5));
  }

  BalancedRibbonFilter(size_t add_count)
      : log2_vshards((uint32_t)FloorLog2((add_count + add_count / 3 + add_count / 5) / (128 * sizeof(CoeffType)))),
        num_slots(GetNumSlots(add_count, log2_vshards)),
        bytes(static_cast<size_t>((num_slots * kFractionalCols + 7) / 8)),
        ptr(new char[bytes]),
        soln(ptr.get(), bytes),
        meta_bytes(BalancedHasher(log2_vshards, nullptr).GetMetadataLength()),
        meta_ptr(new char[meta_bytes]),
        hasher(log2_vshards, meta_ptr.get()) {}

  void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end) {
    BalancedBanding b(log2_vshards);
    b.BalancerAddRange(keys.begin() + start, keys.begin() + end);
    if (!b.Balance(num_slots)) {
      fprintf(stderr, "Failed!\n");
      return;
    }
    soln.BackSubstFrom(b);
    memcpy(meta_ptr.get(), b.GetMetadata(), b.GetMetadataLength());
  }
  bool Contain(uint64_t key) const {
    return soln.FilterQuery(key, hasher);
  }
  size_t SizeInBytes() const {
    return bytes + meta_bytes;
  }
};

template <typename CoeffType, uint32_t kNumColumns, uint32_t kMinPctOverhead, uint32_t kMilliBitsPerKey>
struct FilterAPI<BalancedRibbonFilter<CoeffType, kNumColumns, kMinPctOverhead, kMilliBitsPerKey>> {
  using Table = BalancedRibbonFilter<CoeffType, kNumColumns, kMinPctOverhead, kMilliBitsPerKey>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Contain(key);
  }
};

template <typename CoeffType, uint32_t kNumColumns, uint32_t kMinPctOverhead, bool kUseSmash = false>
class StandardRibbonFilter {
  using TS = RibbonTS<CoeffType, /*kHomog*/ false, kNumColumns, kUseSmash>;
  IMPORT_RIBBON_IMPL_TYPES(TS);

  size_t num_slots;

  size_t bytes;
  unique_ptr<char[]> ptr;
  InterleavedSoln soln;
  Hasher hasher;
public:
  static constexpr double kFractionalCols =
    kNumColumns == 0 ? 7.7 : kNumColumns;

  static double GetNumSlots(size_t add_count) {
    double overhead;
    if (sizeof(CoeffType) == 8) {
      overhead = -0.0251 + std::log(1.0 * add_count) * 1.4427 * 0.0083;
    } else {
      assert(sizeof(CoeffType) == 16);
      overhead = -0.0176 + std::log(1.0 * add_count) * 1.4427 * 0.0038;
    }
    overhead = std::max(overhead, 0.01 * kMinPctOverhead);
    return InterleavedSoln::RoundUpNumSlots((size_t)(add_count + overhead * add_count));
  }

  StandardRibbonFilter(size_t add_count)
      : num_slots(GetNumSlots(add_count)),
        bytes(static_cast<size_t>((num_slots * kFractionalCols + 7) / 8)),
        ptr(new char[bytes]),
        soln(ptr.get(), bytes)
        {}

  void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end) {
    Banding b(num_slots);
    if (!b.AddRange(keys.begin() + start, keys.begin() + end)) {
      fprintf(stderr, "Failed!\n");
      return;
    }
    soln.BackSubstFrom(b);
  }
  bool Contain(uint64_t key) const {
    return soln.FilterQuery(key, hasher);
  }
  size_t SizeInBytes() const {
    return bytes;
  }
};

template <typename CoeffType, uint32_t kNumColumns, uint32_t kMinPctOverhead, bool kUseSmash>
struct FilterAPI<StandardRibbonFilter<CoeffType, kNumColumns, kMinPctOverhead, kUseSmash>> {
  using Table = StandardRibbonFilter<CoeffType, kNumColumns, kMinPctOverhead, kUseSmash>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Contain(key);
  }
};



template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_naive::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_naive::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};


template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_lowmem::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_lowmem::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_naive4wise::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_naive4wise::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};
template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_lowmem4wise::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_lowmem4wise::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};




class MortonFilter {
    Morton3_8* filter;
    size_t size;
public:
    MortonFilter(const size_t size) {
        filter = new Morton3_8((size_t) (size / 0.95) + 64);
        // filter = new Morton3_8((size_t) (2.1 * size) + 64);
        this->size = size;
    }
    ~MortonFilter() {
        delete filter;
    }
    void Add(uint64_t key) {
        filter->insert(key);
    }
    void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end) {
        size_t size = end - start;
        ::std::vector<uint64_t> k(size);
        ::std::vector<bool> status(size);
        for (size_t i = start; i < end; i++) {
            k[i - start] = keys[i];
        }
        // TODO return value and status is ignored currently
        filter->insert_many(k, status, size);
    }
    inline bool Contain(uint64_t &item) {
        return filter->likely_contains(item);
    };
    size_t SizeInBytes() const {
        // according to morton_sample_configs.h:
        // Morton3_8 - 3-slot buckets with 8-bit fingerprints: 11.7 bits/item
        // (load factor = 0.95)
        // so in theory we could just hardcode the size here,
        // and don't measure it
        // return (size_t)((size * 11.7) / 8);

        return filter->SizeInBytes();
    }
};

template<>
struct FilterAPI<MortonFilter> {
    using Table = MortonFilter;
    static Table ConstructFromAddCount(size_t add_count) {
        return Table(add_count);
    }
    static void Add(uint64_t key, Table* table) {
        table->Add(key);
    }
    static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
        table->AddAll(keys, start, end);
    }
    static void Remove(uint64_t, Table *) {
        throw std::runtime_error("Unsupported");
    }
    CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, Table * table) {
        return table->Contain(key);
    }
};

class XorSingle {
public:
    xor8_s filter; // let us expose the struct. to avoid indirection
    explicit XorSingle(const size_t size) {
        if (!xor8_allocate(size, &filter)) {
            throw ::std::runtime_error("Allocation failed");
        }
    }
    ~XorSingle() {
        xor8_free(&filter);
    }
    bool AddAll(const uint64_t* data, const size_t start, const size_t end) {
        return xor8_buffered_populate(data + start, end - start, &filter);
    }
    inline bool Contain(uint64_t &item) const {
        return xor8_contain(item, &filter);
    }
    inline size_t SizeInBytes() const {
        return xor8_size_in_bytes(&filter);
    }
    XorSingle(XorSingle && o) : filter(o.filter)  {
        o.filter.fingerprints = nullptr; // we take ownership for the data
    }
private:
    XorSingle(const XorSingle & o) = delete;
};


template<>
struct FilterAPI<XorSingle> {
    using Table = XorSingle;
    static Table ConstructFromAddCount(size_t add_count) {
        return Table(add_count);
    }
    static void Add(uint64_t, Table*) {
        throw std::runtime_error("Unsupported");
    }
    static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
        table->AddAll(keys.data(), start, end);
    }
    static void Remove(uint64_t, Table *) {
        throw std::runtime_error("Unsupported");
    }
    CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
        // some compilers are not smart enough to do the inlining properly
        return xor8_contain(key, & table->filter);
    }
};

class BinaryFuseSingle {
public:
    binary_fuse8_t filter; // let us expose the struct. to avoid indirection
    explicit BinaryFuseSingle(const size_t size) {
        if (!binary_fuse8_allocate(size, &filter)) {
            throw ::std::runtime_error("Allocation failed");
        }
    }
    ~BinaryFuseSingle() {
        binary_fuse8_free(&filter);
    }
    bool AddAll(const uint64_t* data, const size_t start, const size_t end) {
        return binary_fuse8_populate(data + start, end - start, &filter);
    }
    inline bool Contain(uint64_t &item) const {
        return binary_fuse8_contain(item, &filter);
    }
    inline size_t SizeInBytes() const {
        return binary_fuse8_size_in_bytes(&filter);
    }
    BinaryFuseSingle(BinaryFuseSingle && o) : filter(o.filter)  {
        o.filter.Fingerprints = nullptr; // we take ownership for the data
    }
private:
    BinaryFuseSingle(const BinaryFuseSingle & o) = delete;
};


template<>
struct FilterAPI<BinaryFuseSingle> {
    using Table = BinaryFuseSingle;
    static Table ConstructFromAddCount(size_t add_count) {
        return Table(add_count);
    }
    static void Add(uint64_t , Table* ) {
        throw std::runtime_error("Unsupported");
    }
    static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
        table->AddAll(keys.data(), start, end);
    }
    static void Remove(uint64_t , Table * ) {
        throw std::runtime_error("Unsupported");
    }
    CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
        // some compilers are not smart enough to do the inlining properly
        return binary_fuse8_contain(key, & table->filter);
    }
};


template<size_t blocksize, int k, typename HashFamily>
struct FilterAPI<SimpleBlockFilter<blocksize,k,HashFamily>> {
  using Table = SimpleBlockFilter<blocksize,k,HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(add_count);
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};


template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<XorFilter<ItemType, FingerprintType, HashFamily>> {
  using Table = XorFilter<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};


template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<naive::XorFilter<ItemType, FingerprintType, HashFamily>> {
  using Table = naive::XorFilter<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<prefetch::XorFilter<ItemType, FingerprintType, HashFamily>> {
  using Table = prefetch::XorFilter<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};


template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<XorFilterPlus<ItemType, FingerprintType, HashFamily>> {
  using Table = XorFilterPlus<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<GcsFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = GcsFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

#ifdef __AVX2__
template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<GQFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = GQFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename HashFamily>
struct FilterAPI<VQFilter<ItemType, HashFamily>> {
  using Table = VQFilter<ItemType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  // for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
    // sometimes fails with the original VQF implementation
    // table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};
#endif

template <typename ItemType, size_t bits_per_item, bool branchless, typename HashFamily>
struct FilterAPI<BloomFilter<ItemType, bits_per_item, branchless, HashFamily>> {
  using Table = BloomFilter<ItemType, bits_per_item, branchless, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys.data(), start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, bool branchless, typename HashFamily>
struct FilterAPI<CountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily>> {
  using Table = CountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, bool branchless, typename HashFamily>
struct FilterAPI<SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily>> {
  using Table = SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Contain(key);
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Contain(key);
  }
};


#endif