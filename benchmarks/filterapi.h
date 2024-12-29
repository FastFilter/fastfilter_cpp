#ifndef FILTERAPI_H
#define FILTERAPI_H
#include <climits>
#include <iomanip>
#include <map>
#include <set>
#include <stdexcept>
#include <stdio.h>
#include <vector>

// morton
#include "binaryfusefilter_singleheader.h"
#include "bloom.h"
#include "compressed_cuckoo_filter.h"
#include "counting_bloom.h"
#include "cuckoo_fuse.h"
#include "cuckoofilter.h"
#include "cuckoofilter_stable.h"
#include "gcs.h"
#include "morton_sample_configs.h"
#include "xor_binary_fuse_filter.h"
#include "xorfilter.h"
#include "xorfilter_plus.h"
#include "xorfilter_singleheader.h"
#ifdef __AVX2__
#include "gqf_cpp.h"
#include "simd-block.h"
#include "vqf_cpp.h"
#endif
#define __PF_AVX512__                                                          \
  (__AVX512BW__ & __AVX512VL__ & __AVX512CD__ & __AVX512DQ__)
#if __PF_AVX512__
#include "prefix/min_pd256.hpp"
#include "tc-shortcut/tc-shortcut.hpp"
#endif
#include "ribbon_impl.h"
#include "simd-block-fixed-fpp.h"

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
// compiler. It is unclear whether we want to benchmark the inlining of
// Contains, as it depends very much on how "contains" is used. So it is maybe
// reasonable to benchmark it without inlining.
//
#define CONTAIN_ATTRIBUTES __attribute__((noinline))

template <typename Table> struct FilterAPI {};

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
struct FilterAPI<CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>> {
  using Table = CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t key, Table *table) {
    if (0 != table->Add(key)) {
      throw logic_error("The filter is too small to hold all of the elements");
    }
  }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }
  static void Remove(uint64_t key, Table *table) { table->Delete(key); }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
struct FilterAPI<
    CuckooFilterStable<ItemType, bits_per_item, TableType, HashFamily>> {
  using Table =
      CuckooFilterStable<ItemType, bits_per_item, TableType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t key, Table *table) {
    if (0 != table->Add(key)) {
      throw logic_error("The filter is too small to hold all of the elements");
    }
  }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }
  static void Remove(uint64_t key, Table *table) { table->Delete(key); }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<CuckooFuseFilter<ItemType, FingerprintType>> {
  using Table = CuckooFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t key, Table *table) {
    if (0 != table->Add(key)) {
      throw logic_error("The filter is too small to hold all of the elements");
    }
  }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }
  static void Remove(uint64_t key, Table *table) { table->Delete(key); }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
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
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return table->Find(key);
  }
};

#endif

#ifdef __AVX2__
template <typename HashFamily> struct FilterAPI<SimdBlockFilter<HashFamily>> {
  using Table = SimdBlockFilter<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(ceil(log2(add_count)));
    return ans;
  }
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
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
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
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
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return table->Find(key);
  }
};

#endif
#if __PF_AVX512__
template <typename HashFamily> struct FilterAPI<TC_shortcut<HashFamily>> {
  using Table = TC_shortcut<HashFamily>;

  static Table ConstructFromAddCount(size_t add_count) {
    constexpr float load = .935;
    return Table(add_count, load);
  }
  static void Add(uint64_t key, Table *table) {
    if (!table->insert(key)) {
      std::cout << table->info() << std::endl;
      throw std::logic_error(table->get_name() +
                             " is too small to hold all of the elements");
    }
  }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }

  static bool Add_attempt(uint64_t key, Table *table) {
    if (!table->insert(key)) {
      std::cout << "load when failed: \t" << table->get_effective_load()
                << std::endl;
      std::cout << table->info() << std::endl;
      return false;
    }
    return true;
  }

  static void Remove(uint64_t key, Table *table) { table->remove(key); }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return table->lookup(key);
  }
};

template <typename Table>
inline size_t get_l2_slots(size_t l1_items,
                           const double overflowing_items_ratio,
                           const float loads[2]) {
  const double expected_items_reaching_next_level =
      l1_items * overflowing_items_ratio;
  size_t slots_in_l2 = (expected_items_reaching_next_level / loads[1]);
  return slots_in_l2;
}

template <>
inline size_t get_l2_slots<cuckoofilter::CuckooFilterStable<u64, 12>>(
    size_t l1_items, const double overflowing_items_ratio,
    const float loads[2]) {
  (void)loads;
  (void)overflowing_items_ratio;
  constexpr auto expected_items95 = 0.0586;
  constexpr auto spare_workload = 0.94;
  constexpr auto safety = 1.08;
  constexpr auto factor95 = safety * expected_items95 / spare_workload;
  const double expected_items_reaching_next_level = l1_items * factor95;
  return expected_items_reaching_next_level;
}

template <>
inline size_t get_l2_slots<TC_shortcut<>>(size_t l1_items,
                                          const double overflowing_items_ratio,
                                          const float loads[2]) {
  (void)loads;
  (void)overflowing_items_ratio;
  constexpr auto safety = 1.08;
  constexpr auto expected_items95 = 0.0586;
  constexpr auto spare_workload = 0.935;
  constexpr auto factor95 = safety * expected_items95 / spare_workload;
  const double expected_items_reaching_next_level = l1_items * factor95;
  size_t slots_in_l2 = std::ceil(expected_items_reaching_next_level);
  return slots_in_l2;
}

template <>
inline size_t
get_l2_slots<SimdBlockFilter<>>(size_t l1_items,
                                const double overflowing_items_ratio,
                                const float loads[2]) {
  const double expected_items_reaching_next_level =
      l1_items * overflowing_items_ratio;
  size_t slots_in_l2 = (expected_items_reaching_next_level / loads[1]);
  return slots_in_l2 * 4;
}

template <>
inline size_t
get_l2_slots<SimdBlockFilterFixed<>>(size_t l1_items,
                                     const double overflowing_items_ratio,
                                     const float loads[2]) {
  const double expected_items_reaching_next_level =
      l1_items * overflowing_items_ratio;
  size_t slots_in_l2 = (expected_items_reaching_next_level / loads[1]);
  return slots_in_l2 * 2;
}

template <typename Table,
          typename HashFamily = hashing::TwoIndependentMultiplyShift>
class Prefix_Filter {
  const size_t filter_max_capacity;
  const size_t number_of_pd;
  Table GenSpare;

  hashing::TwoIndependentMultiplyShift Hasher, H0;
  __m256i *pd_array;
  size_t cap[2] = {0};
  static double constexpr overflowing_items_ratio = 0.0586;

public:
  Prefix_Filter(size_t max_items, const float loads[2])
      : filter_max_capacity(max_items),
        number_of_pd(
            std::ceil(1.0 * max_items / (min_pd::MAX_CAP0 * loads[0]))),
        GenSpare(FilterAPI<Table>::ConstructFromAddCount(
            get_l2_slots<Table>(max_items, overflowing_items_ratio, loads))),
        Hasher(), H0() {

    int ok = posix_memalign((void **)&pd_array, 32, 32 * number_of_pd);
    if (ok != 0) {
      std::cout << "Space allocation failed!" << std::endl;
      assert(false);
      exit(-3);
    }

    constexpr uint64_t pd256_plus_init_header =
        (((INT64_C(1) << min_pd::QUOTS) - 1) << 6) | 32;
    for (size_t i = 0; i < number_of_pd; i++) {
      pd_array[i] = __m256i{pd256_plus_init_header, 0, 0, 0};
    }
  }

  ~Prefix_Filter() { free(pd_array); }

  __attribute__((always_inline)) inline static constexpr uint32_t
  reduce32(uint32_t hash, uint32_t n) {
    // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
    return (uint32_t)(((uint64_t)hash * n) >> 32);
  }

  __attribute__((always_inline)) inline static constexpr uint16_t
  fixed_reduce(uint16_t hash) {
    // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
    return (uint16_t)(((uint32_t)hash * 6400) >> 16);
  }

  inline auto Find(const u64 &item) const -> bool {
    const u64 s = H0(item);
    uint32_t out1 = s >> 32u, out2 = s;
    const uint32_t pd_index = reduce32(out1, (uint32_t)number_of_pd);
    const uint16_t qr = fixed_reduce(out2);
    const int64_t quot = qr >> 8;
    const uint8_t rem = qr;
    // return min_pd::pd_find_25(quot, rem, &pd_array[pd_index]);
    // return (!min_pd::cmp_qr1(qr, &pd_array[pd_index])) ?
    // min_pd::pd_find_25(quot, rem, &pd_array[pd_index])
    return (!min_pd::cmp_qr1(qr, &pd_array[pd_index]))
               ? min_pd::find_core(quot, rem, &pd_array[pd_index])
               : incSpare_lookup(pd_index, qr);
  }

  inline auto incSpare_lookup(size_t pd_index, u16 qr) const -> bool {
    const u64 data = (pd_index << 13u) | qr;
    return FilterAPI<Table>::Contain(data, &GenSpare);
  }

  inline void incSpare_add(size_t pd_index, const min_pd::add_res &a_info) {
    cap[1]++;
    u16 qr = (((u16)a_info.quot) << 8u) | a_info.rem;
    const u64 data = (pd_index << 13u) | qr;
    return FilterAPI<Table>::Add(data, &GenSpare);
  }

  void Add(const u64 &item) {
    const u64 s = H0(item);
    constexpr u64 full_mask = (1ULL << 55);
    uint32_t out1 = s >> 32u, out2 = s;

    const uint32_t pd_index = reduce32(out1, (uint32_t)number_of_pd);

    auto pd = pd_array + pd_index;
    const uint64_t header = reinterpret_cast<const u64 *>(pd)[0];
    const bool not_full = !(header & full_mask);

    const uint16_t qr = fixed_reduce(out2);
    const int64_t quot = qr >> 8;
    const uint8_t rem = qr;

    if (not_full) {
      cap[0]++;
      assert(!min_pd::is_pd_full(pd));
      size_t end = min_pd::pd_select64(header >> 6, quot);
      const size_t h_index = end + 6;
      const u64 mask = _bzhi_u64(-1, h_index);
      const u64 lo = header & mask;
      const u64 hi = ((header & ~mask) << 1u); // & h_mask;
      assert(!(lo & hi));
      const u64 h7 = lo | hi;
      memcpy(pd, &h7, 7);

      const size_t body_index = end - quot;
      min_pd::body_add_case0_avx(body_index, rem, pd);
      assert(min_pd::find_core(quot, rem, pd));
      assert(Find(item));
      return;
    } else {
      auto add_res = min_pd::new_pd_swap_short(quot, rem, pd);
      incSpare_add(pd_index, add_res);
      assert(Find(item));
    }
  }

  size_t SizeInBytes() const {
    size_t l1 = sizeof(__m256i) * number_of_pd;
    size_t l2 = GenSpare.SizeInBytes();
    auto res = l1 + l2;
    return res;
  }
};

template <typename filterTable> struct FilterAPI<Prefix_Filter<filterTable>> {
  using Table = Prefix_Filter<filterTable>;

  static Table ConstructFromAddCount(size_t add_count) {
    constexpr float loads[2] = {.95, .95};
    return Table(add_count, loads);
  }

  static void Add(u64 key, Table *table) { table->Add(key); }

  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }

  static void Remove(u64, Table *) { throw std::runtime_error("Unsupported"); }

  CONTAIN_ATTRIBUTES static bool Contain(u64 key, const Table *table) {
    return table->Find(key);
  }
};

#endif

#ifdef __SSE4_1__
template <typename HashFamily>
struct FilterAPI<SimdBlockFilterFixed16<HashFamily>> {
  using Table = SimdBlockFilterFixed16<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(add_count);
    return ans;
  }
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return table->Find(key);
  }
};
#endif

template <typename ItemType, typename FingerprintType>
struct FilterAPI<XorFilter<ItemType, FingerprintType>> {
  using Table = XorFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename CoeffType, bool kHomog, uint32_t kNumColumns,
          bool kSmash = false>
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

  static Hash HashFn(const Hash &input, Seed raw_seed) {
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

template <typename CoeffType, uint32_t kNumColumns,
          uint32_t kMilliBitsPerKey = 7700>
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
    double overhead =
        (4.0 + kFractionalCols * 0.25) / (8.0 * sizeof(CoeffType));
    return 1.0 + overhead;
  }

  HomogRibbonFilter(size_t add_count)
      : num_slots(InterleavedSoln::RoundUpNumSlots(
            (size_t)(GetBestOverheadFactor() * add_count))),
        bytes(static_cast<size_t>((num_slots * kFractionalCols + 7) / 8)),
        ptr(new char[bytes]), soln(ptr.get(), bytes) {}

  void AddAll(const vector<uint64_t> &keys, const size_t start,
              const size_t end) {
    Banding b(num_slots);
    (void)b.AddRange(keys.begin() + start, keys.begin() + end);
    soln.BackSubstFrom(b);
  }
  bool Contain(uint64_t key) const { return soln.FilterQuery(key, hasher); }
  size_t SizeInBytes() const { return bytes; }
};

template <typename CoeffType, uint32_t kNumColumns, uint32_t kMilliBitsPerKey>
struct FilterAPI<HomogRibbonFilter<CoeffType, kNumColumns, kMilliBitsPerKey>> {
  using Table = HomogRibbonFilter<CoeffType, kNumColumns, kMilliBitsPerKey>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return table->Contain(key);
  }
};

template <typename CoeffType, uint32_t kNumColumns, uint32_t kMinPctOverhead,
          uint32_t kMilliBitsPerKey = 7700>
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
    return InterleavedSoln::RoundUpNumSlots(
        (size_t)(add_count + overhead * add_count + add_per_vshard / 5));
  }

  BalancedRibbonFilter(size_t add_count)
      : log2_vshards(
            (uint32_t)FloorLog2((add_count + add_count / 3 + add_count / 5) /
                                (128 * sizeof(CoeffType)))),
        num_slots(GetNumSlots(add_count, log2_vshards)),
        bytes(static_cast<size_t>((num_slots * kFractionalCols + 7) / 8)),
        ptr(new char[bytes]), soln(ptr.get(), bytes),
        meta_bytes(BalancedHasher(log2_vshards, nullptr).GetMetadataLength()),
        meta_ptr(new char[meta_bytes]), hasher(log2_vshards, meta_ptr.get()) {}

  void AddAll(const vector<uint64_t> &keys, const size_t start,
              const size_t end) {
    BalancedBanding b(log2_vshards);
    b.BalancerAddRange(keys.begin() + start, keys.begin() + end);
    if (!b.Balance(num_slots)) {
      fprintf(stderr, "Failed!\n");
      return;
    }
    soln.BackSubstFrom(b);
    memcpy(meta_ptr.get(), b.GetMetadata(), b.GetMetadataLength());
  }
  bool Contain(uint64_t key) const { return soln.FilterQuery(key, hasher); }
  size_t SizeInBytes() const { return bytes + meta_bytes; }
};

template <typename CoeffType, uint32_t kNumColumns, uint32_t kMinPctOverhead,
          uint32_t kMilliBitsPerKey>
struct FilterAPI<BalancedRibbonFilter<CoeffType, kNumColumns, kMinPctOverhead,
                                      kMilliBitsPerKey>> {
  using Table = BalancedRibbonFilter<CoeffType, kNumColumns, kMinPctOverhead,
                                     kMilliBitsPerKey>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return table->Contain(key);
  }
};

template <typename CoeffType, uint32_t kNumColumns, uint32_t kMinPctOverhead,
          bool kUseSmash = false>
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
    return InterleavedSoln::RoundUpNumSlots(
        (size_t)(add_count + overhead * add_count));
  }

  StandardRibbonFilter(size_t add_count)
      : num_slots(GetNumSlots(add_count)),
        bytes(static_cast<size_t>((num_slots * kFractionalCols + 7) / 8)),
        ptr(new char[bytes]), soln(ptr.get(), bytes) {}

  void AddAll(const vector<uint64_t> &keys, const size_t start,
              const size_t end) {
    Banding b(num_slots);
    if (!b.AddRange(keys.begin() + start, keys.begin() + end)) {
      fprintf(stderr, "Failed!\n");
      return;
    }
    soln.BackSubstFrom(b);
  }
  bool Contain(uint64_t key) const { return soln.FilterQuery(key, hasher); }
  size_t SizeInBytes() const { return bytes; }
};

template <typename CoeffType, uint32_t kNumColumns, uint32_t kMinPctOverhead,
          bool kUseSmash>
struct FilterAPI<
    StandardRibbonFilter<CoeffType, kNumColumns, kMinPctOverhead, kUseSmash>> {
  using Table =
      StandardRibbonFilter<CoeffType, kNumColumns, kMinPctOverhead, kUseSmash>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return table->Contain(key);
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<
    xorbinaryfusefilter_naive::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table =
      xorbinaryfusefilter_naive::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_lowmem::XorBinaryFuseFilter<
    ItemType, FingerprintType>> {
  using Table =
      xorbinaryfusefilter_lowmem::XorBinaryFuseFilter<ItemType,
                                                      FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_naive4wise::XorBinaryFuseFilter<
    ItemType, FingerprintType>> {
  using Table =
      xorbinaryfusefilter_naive4wise::XorBinaryFuseFilter<ItemType,
                                                          FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};
template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_lowmem4wise::XorBinaryFuseFilter<
    ItemType, FingerprintType>> {
  using Table =
      xorbinaryfusefilter_lowmem4wise::XorBinaryFuseFilter<ItemType,
                                                           FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

class MortonFilter {
  Morton3_8 *filter;
  size_t size;

public:
  MortonFilter(const size_t size) {
    filter = new Morton3_8((size_t)(size / 0.95) + 64);
    // filter = new Morton3_8((size_t) (2.1 * size) + 64);
    this->size = size;
  }
  ~MortonFilter() { delete filter; }
  void Add(uint64_t key) { filter->insert(key); }
  void AddAll(const vector<uint64_t> &keys, const size_t start,
              const size_t end) {
    size_t size = end - start;
    ::std::vector<uint64_t> k(size);
    ::std::vector<bool> status(size);
    for (size_t i = start; i < end; i++) {
      k[i - start] = keys[i];
    }
    // TODO return value and status is ignored currently
    filter->insert_many(k, status, size);
  }
  inline bool Contain(uint64_t &item) { return filter->likely_contains(item); };
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

template <> struct FilterAPI<MortonFilter> {
  using Table = MortonFilter;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, Table *table) {
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
  ~XorSingle() { xor8_free(&filter); }
  bool AddAll(uint64_t *data, const size_t start, const size_t end) {
    return xor8_buffered_populate(data + start, end - start, &filter);
  }
  inline bool Contain(uint64_t &item) const {
    return xor8_contain(item, &filter);
  }
  inline size_t SizeInBytes() const { return xor8_size_in_bytes(&filter); }
  XorSingle(XorSingle &&o) : filter(o.filter) {
    o.filter.fingerprints = nullptr; // we take ownership for the data
  }

private:
  XorSingle(const XorSingle &o) = delete;
};

template <> struct FilterAPI<XorSingle> {
  using Table = XorSingle;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys.data(), start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    // some compilers are not smart enough to do the inlining properly
    return xor8_contain(key, &table->filter);
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
  ~BinaryFuseSingle() { binary_fuse8_free(&filter); }
  bool AddAll(uint64_t *data, const size_t start, const size_t end) {
    return binary_fuse8_populate(data + start, end - start, &filter);
  }
  inline bool Contain(uint64_t &item) const {
    return binary_fuse8_contain(item, &filter);
  }
  inline size_t SizeInBytes() const {
    return binary_fuse8_size_in_bytes(&filter);
  }
  BinaryFuseSingle(BinaryFuseSingle &&o) : filter(o.filter) {
    o.filter.Fingerprints = nullptr; // we take ownership for the data
  }

private:
  BinaryFuseSingle(const BinaryFuseSingle &o) = delete;
};

template <> struct FilterAPI<BinaryFuseSingle> {
  using Table = BinaryFuseSingle;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys.data(), start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    // some compilers are not smart enough to do the inlining properly
    return binary_fuse8_contain(key, &table->filter);
  }
};

template <size_t blocksize, int k, typename HashFamily>
struct FilterAPI<SimpleBlockFilter<blocksize, k, HashFamily>> {
  using Table = SimpleBlockFilter<blocksize, k, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(add_count);
    return ans;
  }
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return table->Find(key);
  }
};

template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<XorFilter<ItemType, FingerprintType, HashFamily>> {
  using Table = XorFilter<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<naive::XorFilter<ItemType, FingerprintType, HashFamily>> {
  using Table = naive::XorFilter<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<prefetch::XorFilter<ItemType, FingerprintType, HashFamily>> {
  using Table = prefetch::XorFilter<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<XorFilterPlus<ItemType, FingerprintType, HashFamily>> {
  using Table = XorFilterPlus<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<GcsFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = GcsFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

#ifdef __AVX2__
template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<GQFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = GQFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }
  static void Remove(uint64_t key, Table *table) { table->Remove(key); }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename HashFamily>
struct FilterAPI<VQFilter<ItemType, HashFamily>> {
  using Table = VQFilter<ItemType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
    // for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
    // sometimes fails with the original VQF implementation
    // table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};
#endif

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily>
struct FilterAPI<BloomFilter<ItemType, bits_per_item, branchless, HashFamily>> {
  using Table = BloomFilter<ItemType, bits_per_item, branchless, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys.data(), start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily>
struct FilterAPI<
    CountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily>> {
  using Table =
      CountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table *table) { table->Remove(key); }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, bool branchless,
          typename HashFamily>
struct FilterAPI<SuccinctCountingBloomFilter<ItemType, bits_per_item,
                                             branchless, HashFamily>> {
  using Table = SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless,
                                            HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<ItemType> &keys, const size_t start,
                     const size_t end, Table *table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table *table) { table->Remove(key); }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<
    SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily>> {
  using Table =
      SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }
  static void Remove(uint64_t key, Table *table) { table->Remove(key); }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return table->Contain(key);
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item,
                                                        HashFamily>> {
  using Table = SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item,
                                                       HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    return Table(add_count);
  }
  static void Add(uint64_t key, Table *table) { table->Add(key); }
  static void AddAll(const vector<uint64_t> &keys, const size_t start,
                     const size_t end, Table *table) {
    for (size_t i = start; i < end; i++) {
      Add(keys[i], table);
    }
  }
  static void Remove(uint64_t key, Table *table) { table->Remove(key); }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table *table) {
    return table->Contain(key);
  }
};

#endif
