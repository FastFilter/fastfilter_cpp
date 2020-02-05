// This benchmark reports on the bulk insert and bulk query rates. It is invoked as:
//
//     ./bulk-insert-and-query.exe 158000
//
// That invocation will test each probabilistic membership container type with 158000
// randomly generated items. It tests bulk Add() from empty to full and Contain() on
// filters with varying rates of expected success. For instance, at 75%, three out of
// every four values passed to Contain() were earlier Add()ed.
//
// Example usage:
//
// for alg in `seq 0 1 14`; do for num in `seq 10 10 200`; do ./bulk-insert-and-query.exe ${num}000000 ${alg}; done; done > results.txt

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
#include "xorfilter.h"
#include "xorfilter_10bit.h"
#include "xorfilter_13bit.h"
#include "xorfilter_10_666bit.h"
#include "xorfilter_2.h"
#include "xorfilter_2n.h"
#include "xorfilter_plus.h"
#include "xorfilter_singleheader.h"
#include "xor_fuse_filter.h"
#include "bloom.h"
#include "counting_bloom.h"
#include "gcs.h"
#ifdef __AVX2__
#include "gqf_cpp.h"
#include "simd-block.h"
#endif
#include "random.h"
#include "simd-block-fixed-fpp.h"
#include "timing.h"
#ifdef __linux__
#include "linux-perf-events.h"
#endif

using namespace std;
using namespace hashing;
using namespace cuckoofilter;
using namespace xorfilter;
using namespace xorfilter2;
using namespace xorfilter2n;
using namespace xorfilter_plus;
using namespace xorfusefilter;
using namespace bloomfilter;
using namespace counting_bloomfilter;
using namespace gcsfilter;
using namespace CompressedCuckoo; // Morton filter namespace
#ifdef __AVX2__
using namespace gqfilter;
#endif

// The number of items sampled when determining the lookup performance
const size_t MAX_SAMPLE_SIZE = 10 * 1000 * 1000;

// The statistics gathered for each table type:
struct Statistics {
  size_t add_count;
  double nanos_per_add;
  double nanos_per_remove;
  // key: percent of queries that were expected to be positive
  map<int, double> nanos_per_finds;
  double false_positive_probabilty;
  double bits_per_item;
};

// Inlining the "contains" which are executed within a tight loop can be both
// very detrimental or very beneficial, and which ways it goes depends on the
// compiler. It is unclear whether we want to benchmark the inlining of Contains,
// as it depends very much on how "contains" is used. So it is maybe reasonable
// to benchmark it without inlining.
//
#define CONTAIN_ATTRIBUTES  __attribute__ ((noinline))

// Output for the first row of the table of results. type_width is the maximum number of
// characters of the description of any table type, and find_percent_count is the number
// of different lookup statistics gathered for each table. This function assumes the
// lookup expected positive probabiilties are evenly distributed, with the first being 0%
// and the last 100%.
string StatisticsTableHeader(int type_width, int find_percent_count) {
  ostringstream os;

  os << string(type_width, ' ');
  os << setw(8) << right << "";
  os << setw(8) << right << "";
  for (int i = 0; i < find_percent_count; ++i) {
    os << setw(8) << "find";
  }
  os << setw(9) << "" << setw(11) << "" << setw(11)
     << "optimal" << setw(8) << "wasted" << setw(8) << "million" << endl;

  os << string(type_width, ' ');
  os << setw(8) << right << "add";
  os << setw(8) << right << "remove";
  for (int i = 0; i < find_percent_count; ++i) {
    os << setw(7)
       << static_cast<int>(100 * i / static_cast<double>(find_percent_count - 1)) << '%';
  }
  os << setw(10) << "ε" << setw(11) << "bits/item" << setw(11)
     << "bits/item" << setw(8) << "space" << setw(8) << "keys";
  return os.str();
}

// Overloading the usual operator<< as used in "std::cout << foo", but for Statistics
template <class CharT, class Traits>
basic_ostream<CharT, Traits>& operator<<(
    basic_ostream<CharT, Traits>& os, const Statistics& stats) {
  os << fixed << setprecision(2) << setw(8) << right
     << stats.nanos_per_add;
  os << fixed << setprecision(2) << setw(8) << right
     << stats.nanos_per_remove;
  for (const auto& fps : stats.nanos_per_finds) {
    os << setw(8) << fps.second;
  }
  // we get some nonsensical result for very small fpps
  if(stats.false_positive_probabilty > 0.0000001) {
    const auto minbits = log2(1 / stats.false_positive_probabilty);
    os << setw(8) << setprecision(4) << stats.false_positive_probabilty * 100 << '%'
       << setw(11) << setprecision(2) << stats.bits_per_item << setw(11) << minbits
       << setw(7) << setprecision(1) << 100 * (stats.bits_per_item / minbits - 1) << '%'
       << setw(8) << setprecision(1) << (stats.add_count / 1000000.);
  } else {
    os << setw(8) << setprecision(4) << stats.false_positive_probabilty * 100 << '%'
       << setw(11) << setprecision(2) << stats.bits_per_item << setw(11) << 64
       << setw(7) << setprecision(1) << 0 << '%'
       << setw(8) << setprecision(1) << (stats.add_count / 1000000.);
  }
  return os;
}

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
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    throw std::runtime_error("Unsupported");
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
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    throw std::runtime_error("Unsupported");
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
    Table ans(ceil(add_count * 8.0 / CHAR_BIT));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t> keys, const size_t start, const size_t end, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void Remove(uint64_t key, Table * table) {
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
    Table ans(ceil(log2(add_count * 8.0 / CHAR_BIT)));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t> keys, const size_t start, const size_t end, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void Remove(uint64_t key, Table * table) {
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
    Table ans(ceil(add_count * 8.0 / CHAR_BIT));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t> keys, const size_t start, const size_t end, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void Remove(uint64_t key, Table * table) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};


template <typename HashFamily>
struct FilterAPI<SimdBlockFilterFixed16<HashFamily>> {
  using Table = SimdBlockFilterFixed16<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(ceil(add_count * 8.0 / CHAR_BIT));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t> keys, const size_t start, const size_t end, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void Remove(uint64_t key, Table * table) {
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
    Table ans(ceil(add_count * 8.0 / CHAR_BIT));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
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
  static void Add(uint64_t key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<XorFuseFilter<ItemType, FingerprintType>> {
  using Table = XorFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
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
    void AddAll(const vector<uint64_t> keys, const size_t start, const size_t end) {
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
    static void AddAll(const vector<uint64_t> keys, const size_t start, const size_t end, Table* table) {
        table->AddAll(keys, start, end);
    }
    static void Remove(uint64_t key, Table * table) {
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
    static void Add(uint64_t key, Table* table) {
        throw std::runtime_error("Unsupported");
    }
    static void AddAll(const vector<uint64_t> keys, const size_t start, const size_t end, Table* table) {
        table->AddAll(keys.data(), start, end);
    }
    static void Remove(uint64_t key, Table * table) {
        throw std::runtime_error("Unsupported");
    }
    CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
        // some compilers are not smart enough to do the inlining properly
        return xor8_contain(key, & table->filter);
    }
};

template<size_t blocksize, int k, typename HashFamily>
struct FilterAPI<SimpleBlockFilter<blocksize,k,HashFamily>> {
  using Table = SimpleBlockFilter<blocksize,k,HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(ceil(add_count * 8.0 / CHAR_BIT));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t> keys, const size_t start, const size_t end, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void Remove(uint64_t key, Table * table) {
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
  static void Add(uint64_t key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename FingerprintStorageType, typename HashFamily>
struct FilterAPI<XorFilter2<ItemType, FingerprintType, FingerprintStorageType, HashFamily>> {
  using Table = XorFilter2<ItemType, FingerprintType, FingerprintStorageType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename HashFamily>
struct FilterAPI<XorFilter10<ItemType, HashFamily>> {
  using Table = XorFilter10<ItemType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename HashFamily>
struct FilterAPI<XorFilter13<ItemType, HashFamily>> {
  using Table = XorFilter13<ItemType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename HashFamily>
struct FilterAPI<XorFilter10_666<ItemType, HashFamily>> {
  using Table = XorFilter10_666<ItemType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename FingerprintStorageType, typename HashFamily>
struct FilterAPI<XorFilter2n<ItemType, FingerprintType, FingerprintStorageType, HashFamily>> {
  using Table = XorFilter2n<ItemType, FingerprintType, FingerprintStorageType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
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
  static void Add(uint64_t key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
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
  static void Add(uint64_t key, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
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
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    throw std::runtime_error("Unsupported");
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
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
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
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
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
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
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
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
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    throw std::runtime_error("Unsupported");
    // table->AddAll(keys, start, end);
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
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    throw std::runtime_error("Unsupported");
    // table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Contain(key);
  }
};

// assuming that first1,last1 and first2, last2 are sorted,
// this tries to find out how many of first1,last1 can be
// found in first2, last2, this includes duplicates
template<class InputIt1, class InputIt2>
size_t match_size_iter(InputIt1 first1, InputIt1 last1,
                          InputIt2 first2, InputIt2 last2) {
    size_t answer = 0;
    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2) {
            ++first1;
        } else  if (*first2 < *first1) {
            ++first2;
        } else {
            answer ++;
            ++first1;
        }
    }
    return answer;
}

template<class InputIt>
size_t count_distinct(InputIt first, InputIt last) {
    if(last  == first) return 0;
    size_t answer = 1;
    auto val = *first;
    first++;

    while (first != last) {
      if(val != *first) ++answer;
      first++;
    }
    return answer;
}

size_t match_size(vector<uint64_t> a,  vector<uint64_t> b, size_t * distincta, size_t * distinctb) {
  // could obviously be accelerated with a Bloom filter
  // But this is surprisingly fast!
  vector<uint64_t> result;
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());
  if(distincta != NULL) *distincta  = count_distinct(a.begin(), a.end());
  if(distinctb != NULL) *distinctb  = count_distinct(b.begin(), b.end());
  return match_size_iter(a.begin(), a.end(),b.begin(), b.end());
}

bool has_duplicates(vector<uint64_t> a) {
  std::sort(a.begin(), a.end());
  return count_distinct(a.begin(), a.end()) < a.size();
}
struct samples {
  double found_probability;
  std::vector<uint64_t> to_lookup_mixed;
  size_t true_match;
  size_t actual_sample_size;
};

typedef struct samples samples_t;

template <typename Table>
Statistics FilterBenchmark(
    size_t add_count, const vector<uint64_t>& to_add, size_t distinct_add, const vector<uint64_t>& to_lookup, size_t distinct_lookup,
    size_t intersectionsize, bool hasduplicates,
    const std::vector<samples_t> & mixed_sets, int seed, bool batchedadd = false, bool remove = false) {
  if (add_count > to_add.size()) {
    throw out_of_range("to_add must contain at least add_count values");
  }


  Table filter = FilterAPI<Table>::ConstructFromAddCount(add_count);
  Statistics result;
#ifdef __linux__
  vector<int> evts;
  evts.push_back(PERF_COUNT_HW_CPU_CYCLES);
  evts.push_back(PERF_COUNT_HW_INSTRUCTIONS);
  evts.push_back(PERF_COUNT_HW_CACHE_MISSES);
  evts.push_back(PERF_COUNT_HW_BRANCH_MISSES);
  LinuxEvents<PERF_TYPE_HARDWARE> unified(evts);
  vector<unsigned long long> results;
  results.resize(evts.size());
  cout << endl;
  unified.start();
#else
   std::cout << "-" << std::flush;
#endif

  // Add values until failure or until we run out of values to add:
  if(batchedadd) {
    std::cout << "batched add" << std::flush;
  } else {
    std::cout << "1-by-1 add" << std::flush;
  }
  auto start_time = NowNanos();
  if(batchedadd) {
    FilterAPI<Table>::AddAll(to_add, 0, add_count, &filter);
  } else {
    for (size_t added = 0; added < add_count; ++added) {
      FilterAPI<Table>::Add(to_add[added], &filter);
    }
  }
  auto time = NowNanos() - start_time;
  std::cout << "\r             \r" << std::flush;
#ifdef __linux__
  unified.end(results);
  printf("add    ");
  printf("cycles: %5.1f/key, instructions: (%5.1f/key, %4.2f/cycle) cache misses: %5.2f/key branch misses: %4.2f/key\n",
    results[0]*1.0/add_count,
    results[1]*1.0/add_count ,
    results[1]*1.0/results[0],
    results[2]*1.0/add_count,
    results[3]*1.0/add_count);
#else
  std::cout << "." << std::flush;
#endif

  // sanity check:
  for (size_t added = 0; added < add_count; ++added) {
    assert(FilterAPI<Table>::Contain(to_add[added], &filter) == 1);
  }

  result.add_count = add_count;
  result.nanos_per_add = static_cast<double>(time) / add_count;
  result.bits_per_item = static_cast<double>(CHAR_BIT * filter.SizeInBytes()) / add_count;
  size_t found_count = 0;

  for (auto t :  mixed_sets) {
    const double found_probability = t.found_probability;
    const auto to_lookup_mixed =  t.to_lookup_mixed ;
    size_t true_match = t.true_match ;

#ifdef __linux__
    unified.start();
#else
    std::cout << "-" << std::flush;
#endif
    const auto start_time = NowNanos();
    found_count = 0;
    for (const auto v : to_lookup_mixed) {
      found_count += FilterAPI<Table>::Contain(v, &filter);
    }
    const auto lookup_time = NowNanos() - start_time;
#ifdef __linux__
    unified.end(results);
    printf("%3.2f%%  ",found_probability);
    printf("cycles: %5.1f/key, instructions: (%5.1f/key, %4.2f/cycle) cache misses: %5.2f/key branch misses: %4.2f/key\n",
      results[0]*1.0/to_lookup_mixed.size(),
      results[1]*1.0/to_lookup_mixed.size(),
      results[1]*1.0/results[0],
      results[2]*1.0/to_lookup_mixed.size(),
      results[3] * 1.0/to_lookup_mixed.size());
#else
    std::cout << "." << std::flush;
#endif

    if (found_count < true_match) {
           cerr << "ERROR: Expected to find at least " << true_match << " found " << found_count << endl;
           cerr << "ERROR: This is a potential bug!" << endl;
    }
    result.nanos_per_finds[100 * found_probability] =
        static_cast<double>(lookup_time) / t.actual_sample_size;
    if (0.0 == found_probability) {
      ////////////////////////////
      // This is obviously technically wrong!!! The assumption is that there is no overlap between the random
      // queries and the random content. This is likely true if your 64-bit values were generated randomly,
      // but not true in general.
      ///////////////////////////
      // result.false_positive_probabilty =
      //    found_count / static_cast<double>(to_lookup_mixed.size());
      if(t.to_lookup_mixed.size() == intersectionsize) {
        cerr << "WARNING: fpp is probably meaningless! " << endl;
      }
      result.false_positive_probabilty = (found_count  - intersectionsize) / static_cast<double>(to_lookup_mixed.size() - intersectionsize);
    }
  }

  // Remove
  result.nanos_per_remove = 0;
  if (remove) {
    std::cout << "1-by-1 remove" << std::flush;
#ifdef __linux__
    unified.start();
#else
    std::cout << "-" << std::flush;
#endif
    start_time = NowNanos();
    for (size_t added = 0; added < add_count; ++added) {
      FilterAPI<Table>::Remove(to_add[added], &filter);
    }
    time = NowNanos() - start_time;
    result.nanos_per_remove = static_cast<double>(time) / add_count;
#ifdef __linux__
    unified.end(results);
    printf("remove ");
    printf("cycles: %5.1f/key, instructions: (%5.1f/key, %4.2f/cycle) cache misses: %5.2f/key branch misses: %4.2f/key\n",
      results[0]*1.0/add_count,
      results[1]*1.0/add_count ,
      results[1]*1.0/results[0],
      results[2]*1.0/add_count,
      results[3]*1.0/add_count);
#else
    std::cout << "." << std::flush;
#endif
  }

#ifndef __linux__
  std::cout << "\r             \r" << std::flush;
#endif

  return result;
}

uint64_t reverseBitsSlow(uint64_t v) {
    // r will be reversed bits of v; first get LSB of v
    uint64_t r = v & 1;
    int s = sizeof(v) * CHAR_BIT - 1; // extra shift needed at end
    for (v >>= 1; v; v >>= 1) {
        r <<= 1;
        r |= v & 1;
        s--;
    }
    r <<= s; // shift when v's highest bits are zero
    return r;
}

void parse_comma_separated(char * c, std::set<int> & answer ) {
    std::stringstream ss(c);
    int i;
    while (ss >> i) {
        answer.insert(i);
        if (ss.peek() == ',')
            ss.ignore();
    }
}


int main(int argc, char * argv[]) {
  std::map<int,std::string> names = {
    // Xor
    {0, "Xor8"}, {1, "Xor12"}, {2, "Xor16"},
    {3, "Xor+8"}, {4, "Xor+16"},
    {5, "Xor10"}, {6, "Xor10.666"},
    {7, "Xor10 (NBitArray)"}, {8, "Xor14 (NBitArray)"}, {9, "Xor8-2^n"},
    // Cuckooo
    {10,"Cuckoo8"}, {11,"Cuckoo12"}, {12,"Cuckoo16"},
    {13,"CuckooSemiSort13"},
    {14, "Cuckoo8-2^n"}, {15, "Cuckoo12-2^n"}, {16, "Cuckoo16-2^n"},
    {17, "CuckooSemiSort13-2^n"},
    // GCS
    {20,"GCS"},
#ifdef __AVX2__
    // CQF
    {30,"CQF"},
#endif
    // Bloom
    {40, "Bloom8"}, {41, "Bloom12" }, {42, "Bloom16"},
    {43, "Bloom8 (addall)"}, {44, "Bloom12 (addall)"}, {45, "Bloom16 (addall)"},
    {46, "BranchlessBloom8 (addall)"},
    {47, "BranchlessBloom12 (addall)"},
    {48, "BranchlessBloom16 (addall)"},
    // Blocked Bloom
    {50, "SimpleBlockedBloom"},
#ifdef __aarch64__
    {51, "BlockedBloom"},
    {52, "BlockedBloom (addall)"},
#elif defined( __AVX2__)
    {51, "BlockedBloom"},
    {52, "BlockedBloom (addall)"},
    {53, "BlockedBloom64"},
#endif
#ifdef __SSE41__
    {54, "BlockedBloom16"},
#endif

    // Counting Bloom
    {60, "CountingBloom10 (addall)"},
    {61, "SuccCountingBloom10 (addall)"},
    {62, "SuccCountBlockBloom10"},
    {63, "SuccCountBlockBloomRank10"},

    {70, "Xor8-singleheader"},
    {80, "Morton"},

    {90, "XorFuse8"},

    // Sort
    {100, "Sort"},
  };

  // Parameter Parsing ----------------------------------------------------------

  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <numberOfEntries> [<algorithmId> [<seed>]]" << endl;
    cout << " numberOfEntries: number of keys, we recommend at least 100000000" << endl;
    cout << " algorithmId: -1 for all default algos, or 0..n to only run this algorithm" << endl;
    cout << " algorithmId: can also be a comma-separated list of non-negative integers" << endl;
    for(auto i : names) {
      cout << "           "<< i.first << " : " << i.second << endl;
    }
    cout << " algorithmId: can also be set to the string 'all' if you want to run them all, including some that are excluded by default" << endl;
    cout << " seed: seed for the PRNG; -1 for random seed (default)" << endl;
    return 1;
  }
  stringstream input_string(argv[1]);
  size_t add_count;
  input_string >> add_count;
  if (input_string.fail()) {
    cerr << "Invalid number: " << argv[1];
    return 2;
  }
  int algorithmId = -1; // -1 is just the default
  std::set<int> algos;
  if (argc > 2) {
      if(strcmp(argv[2],"all") == 0) {
         for(auto i : names) {// we add all the named algos.
           algos.insert(i.first);
         }
      } else if(strstr(argv[2],",") != NULL) {
        // we have a list of algos
        algorithmId = 9999999; // disabling
        parse_comma_separated(argv[2], algos);
        if(algos.size() == 0) {
           cerr<< " no algo selected " << endl;
           return -3;
        }
      } else {
        // we select just one
        stringstream input_string_2(argv[2]);
        input_string_2 >> algorithmId;
        if (input_string_2.fail()) {
            cerr << "Invalid number: " << argv[2];
            return 2;
        }
      }
  }
  int seed = -1;
  if (argc > 3) {
      stringstream input_string_3(argv[3]);
      input_string_3 >> seed;
      if (input_string_3.fail()) {
          cerr << "Invalid number: " << argv[3];
          return 2;
      }
  }
  size_t actual_sample_size = MAX_SAMPLE_SIZE;
  if (actual_sample_size > add_count) {
    actual_sample_size = add_count;
  }

  // Generating Samples ----------------------------------------------------------

  vector<uint64_t> to_add = seed == -1 ?
      GenerateRandom64Fast(add_count, rand()) :
      GenerateRandom64Fast(add_count, seed);
  vector<uint64_t> to_lookup = seed == -1 ?
      GenerateRandom64Fast(actual_sample_size, rand()) :
      GenerateRandom64Fast(actual_sample_size, seed + add_count);

  if (seed >= 0 && seed < 64) {
    // 0-64 are special seeds
    uint rotate = seed;
    cout << "Using sequential ordering rotated by " << rotate << endl;
    for(uint64_t i = 0; i < to_add.size(); i++) {
        to_add[i] = xorfilter::rotl64(i, rotate);
    }
    for(uint64_t i = 0; i < to_lookup.size(); i++) {
        to_lookup[i] = xorfilter::rotl64(i + to_add.size(), rotate);
    }
  } else if (seed >= 64 && seed < 128) {
    // 64-127 are special seeds
    uint rotate = seed - 64;
    cout << "Using sequential ordering rotated by " << rotate << " and reversed bits " << endl;
    for(uint64_t i = 0; i < to_add.size(); i++) {
        to_add[i] = reverseBitsSlow(xorfilter::rotl64(i, rotate));
    }
    for(uint64_t i = 0; i < to_lookup.size(); i++) {
        to_lookup[i] = reverseBitsSlow(xorfilter::rotl64(i + to_add.size(), rotate));
    }
  }

  assert(to_lookup.size() == actual_sample_size);
  size_t distinct_lookup;
  size_t distinct_add;
  std::cout << "checking match size... " << std::flush;
  size_t intersectionsize = match_size(to_lookup, to_add, &distinct_lookup, & distinct_add);
  std::cout << "\r                       \r" << std::flush;

  bool hasduplicates = false;
  if(intersectionsize > 0) {
    cout << "WARNING: Out of the lookup table, "<< intersectionsize<< " ("<<intersectionsize * 100.0 / to_lookup.size() << "%) of values are present in the filter." << endl;
    hasduplicates = true;
  }

  if(distinct_lookup != to_lookup.size()) {
    cout << "WARNING: Lookup contains "<< (to_lookup.size() - distinct_lookup)<<" duplicates." << endl;
    hasduplicates = true;
  }
  if(distinct_add != to_add.size()) {
    cout << "WARNING: Filter contains "<< (to_add.size() - distinct_add) << " duplicates." << endl;
    hasduplicates = true;
  }

  if (actual_sample_size > to_lookup.size()) {
    std::cerr << "actual_sample_size = "<< actual_sample_size << std::endl;
    throw out_of_range("to_lookup must contain at least actual_sample_size values");
  }

  std::vector<samples_t> mixed_sets;

  for (const double found_probability : {0.0, 0.25, 0.50, 0.75, 1.00}) {
    std::cout << "generating samples with probability " << found_probability <<" ... " << std::flush;

    struct samples thisone;
    thisone.found_probability = found_probability;
    thisone.actual_sample_size = actual_sample_size;
    uint64_t mixingseed = seed == -1 ? random() : seed;
    thisone.to_lookup_mixed = DuplicateFreeMixIn(&to_lookup[0], &to_lookup[actual_sample_size], &to_add[0],
    &to_add[add_count], found_probability, mixingseed);
    assert(thisone.to_lookup_mixed.size() == actual_sample_size);
    thisone.true_match = match_size(thisone.to_lookup_mixed,to_add, NULL, NULL);
    double trueproba =  thisone.true_match /  static_cast<double>(actual_sample_size) ;
    double bestpossiblematch = fabs(round(found_probability * actual_sample_size) / static_cast<double>(actual_sample_size) - found_probability);
    double tolerance = bestpossiblematch > 0.01 ? bestpossiblematch : 0.01;
    double probadiff = fabs(trueproba - found_probability);
    if(probadiff >= tolerance) {
      cerr << "WARNING: You claim to have a find proba. of " << found_probability << " but actual is " << trueproba << endl;
      return EXIT_FAILURE;
    }
    mixed_sets.push_back(thisone);
    std::cout << "\r                                                                                         \r"  << std::flush;
  }
  constexpr int NAME_WIDTH = 32;
  cout << StatisticsTableHeader(NAME_WIDTH, 5) << endl;

  // Algorithms ----------------------------------------------------------
  int a;

  // Xor ----------------------------------------------------------
  a = 0;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter<uint64_t, uint8_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 1;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter2<uint64_t, uint32_t, UInt12Array, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 2;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter<uint64_t, uint16_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilterPlus<uint64_t, uint8_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 4;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilterPlus<uint64_t, uint16_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 5;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter10<uint64_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 6;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter10_666<uint64_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 7;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter2<uint64_t, uint16_t, NBitArray<uint16_t, 10>, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 8;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter2<uint64_t, uint16_t, NBitArray<uint16_t, 14>, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 9;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter2n<uint64_t, uint8_t, UIntArray<uint8_t>, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  // Cuckoo ----------------------------------------------------------
  a = 10;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 8, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 11;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 12, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 12;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 16, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 13;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 13, PackedTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 14;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 8, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 15;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 12, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 16;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 16, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 17;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 13, PackedTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  // GCS ----------------------------------------------------------
  a = 20;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          GcsFilter<uint64_t, 8, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  // CQF ----------------------------------------------------------
#ifdef __AVX2__
  a = 30;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          GQFilter<uint64_t, 8, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
#endif

  // Bloom ----------------------------------------------------------
  a = 40;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 8, false, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 41;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 12, false, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 42;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 16, false, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 43;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 8, false, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 44;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 12, false, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 45;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 16, false, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 46;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 8, true, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 47;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 12, true, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 48;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 16, true, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  a = 48;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 16, true, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  // Blocked Bloom ----------------------------------------------------------
  a = 50;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          SimpleBlockFilter<8, 8, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
#ifdef __aarch64__
  a = 51;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed<SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 52;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed<SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
#endif
#ifdef __AVX2__
  a = 51;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed<SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 52;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed<SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 53;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
        auto cf = FilterBenchmark<SimdBlockFilterFixed64<SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
#endif
#ifdef __SSE41__
  a = 54;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed16<SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
#endif

  // Counting Bloom ----------------------------------------------------------
  a = 60;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CountingBloomFilter<uint64_t, 10, true, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 61;
  if (algorithmId == a  || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          SuccinctCountingBloomFilter<uint64_t, 10, true, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 62;
  if (algorithmId == a  || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          SuccinctCountingBlockedBloomFilter<uint64_t, 10, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 63;
  if (algorithmId == a  || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          SuccinctCountingBlockedBloomRankFilter<uint64_t, 10, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  a = 70;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorSingle>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  a = 80;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          MortonFilter>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  // Xor Fuse Filter ----------------------------------------------------------
  a = 90;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFuseFilter<uint64_t, uint8_t>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  // Sort ----------------------------------------------------------
  a = 100;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto start_time = NowNanos();
      std::sort(to_add.begin(), to_add.end());
      const auto sort_time = NowNanos() - start_time;
      std::cout << "Sort time: " << sort_time / to_add.size() << " ns/key\n";
  }


}
