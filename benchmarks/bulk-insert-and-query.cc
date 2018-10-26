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

#include "cuckoofilter.h"
#include "cuckoofilter_stable.h"
#include "xorfilter.h"
#include "xorfilter_10bit.h"
#include "xorfilter_13bit.h"
#include "xorfilter_10_666bit.h"
#include "xorfilter_2.h"
#include "xorfilter_2n.h"
#include "xorfilter_plus.h"
#include "bloom.h"
#include "gcs.h"
#ifdef __AVX2__
#include "gqf_cpp.h"
#endif
#include "random.h"
#ifdef __AVX2__
#include "simd-block.h"
#include "simd-block-fixed-fpp.h"
#endif
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
using namespace bloomfilter;
using namespace gcsfilter;
#ifdef __AVX2__
using namespace gqfilter;
#endif

// The number of items sampled when determining the lookup performance
const size_t SAMPLE_SIZE = 10 * 1000 * 1000;

// The statistics gathered for each table type:
struct Statistics {
  size_t add_count;
  double nanos_per_add;
  map<int, double> nanos_per_finds; // The key is the percent of queries that were expected
                                   // to be positive
  double false_positive_probabilty;
  double bits_per_item;
};

//
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
  os << setw(12) << right << "";
  for (int i = 0; i < find_percent_count; ++i) {
    os << setw(8) << "find";
  }
  os << setw(9) << "" << setw(11) << "" << setw(11)
     << "optimal" << setw(8) << "wasted" << setw(8) << "million" << endl;

  os << string(type_width, ' ');
  os << setw(12) << right << "ns/add";
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
  os << fixed << setprecision(2) << setw(12) << right
     << stats.nanos_per_add;
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
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
  static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
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

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};
#endif

template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<BloomFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = BloomFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }

  CONTAIN_ATTRIBUTES
  static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
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
    const std::vector<samples_t> & mixed_sets, int seed, bool batchedadd = false) {
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
  auto start_time = NowNanos();
  if(batchedadd) {
    std::cout << "batched add" << std::flush;
    // for the XorFilter
    FilterAPI<Table>::AddAll(to_add, 0, add_count, &filter);
  } else {
    std::cout << "1-by-1 add" << std::flush;
    for (size_t added = 0; added < add_count; ++added) {
      FilterAPI<Table>::Add(to_add[added], &filter);
    }
  }
  std::cout << "\r             \r" << std::flush;

  // sanity check:
  for (size_t added = 0; added < add_count; ++added) {
    assert(FilterAPI<Table>::Contain(to_add[added], &filter) == 1);
  }
  auto time = NowNanos() - start_time;

#ifdef __linux__
    unified.end(results);
    printf("cycles: %10.zu (%10.3f per key) instructions: %10.zu (%10.3f per key, %10.3f per cycle) cache misses: %10.zu (%10.3f per key) branch misses: %10.zu (%10.3f per key)\n",
      (size_t)results[0], results[0]*1.0/add_count, (size_t)results[1], results[1]*1.0/add_count , results[1]*1.0/results[0], (size_t)results[2], results[2]*1.0/add_count,
      (size_t)results[3], results[3] * 1.0/add_count);
#else
   std::cout << "." << std::flush;
#endif

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
    printf("cycles: %10.zu (%10.3f per key) instructions: %10.zu (%10.3f per key, %10.3f per cycle) cache misses: %10.zu (%10.3f per key) branch misses: %10.zu (%10.3f per key)\n",
      (size_t)results[0], results[0]*1.0/to_lookup_mixed.size(), (size_t)results[1], results[1]*1.0/to_lookup_mixed.size() , results[1]*1.0/results[0], (size_t)results[2], results[2]*1.0/to_lookup_mixed.size(),
      (size_t)results[3], results[3] * 1.0/to_lookup_mixed.size());
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
#ifndef __linux__
  std::cout <<"\r"; // roll back the line
#endif
  return result;
}

void fixEndian(uint64_t* longArray, uint64_t byteCount) {
    uint8_t* byteArray = (uint8_t*) longArray;
    uint64_t l=0;
    uint64_t b=0;
    while(b<byteCount) {
        uint64_t x = 0;
        for(int i=0; i<8; i++) {
            x = (x << 8) | byteArray[b++];
        }
        longArray[l++] = x;
    }
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

/*
#define MUL 1625L
#define MUL2 (MUL*MUL)
// (1<<64) / MUL
#define INVMUL 11351842506898185L
// (1<<64) / MUL2
#define INVMUL2 6985749235014L
int main() {
    printf("start\n");
    for (int a = 0; a < MUL; a++) {
        for (int b = 0; b < MUL; b++) {
            for (int c = 0; c < MUL; c++) {
                uint32_t x = a * MUL2 + b * MUL + c;
                int aa = (int) (((__uint128_t) x * (INVMUL2 + 1)) >> 64);
                if (aa != a) {
                    printf("wrong a");
                    return -1;
                }
                int bb = (int) (((__uint128_t) x * (INVMUL + 1)) >> 64);
                int rb = bb % MUL;
                if (rb != b) {
                    printf("wrong b");
                    return -1;
                }
                int expected = (a + b + c) % MUL;
                int got = (aa + bb + x) % MUL;
                if (expected != got) {
                    printf("wrong modulo");
                    return -1;
                }
            }
        }
    }
    printf("end\n");
    return 0;
}
*/

int main(int argc, char * argv[]) {
  std::map<int,std::string> names = {{0,"Xor8"},{1,"Xor12"},
   {2,"Xor16"}, {3,"Cuckoo8"}, {4,"Cuckoo12"},
   {5,"Cuckoo16"}, {6,"CuckooSemiSort13" }, {7,"Bloom8"},
   {8,"Bloom12" }, {9,"Bloom16"}, {10,"BlockedBloom"},
   {11,"sort"}, {12,"Xor+8"}, {13,"Xor+16"},
   {14,"GCS"}, {15,"CQF"}, {25, "Xor10"}, {37,"Bloom8 (addall)"},
   {38,"Bloom12 (addall)"},{39,"Bloom16 (addall)"},
   {40,"BlockedBloom (addall)"}
  };
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <numberOfEntries> [<algorithmId> [<seed>]]" << endl;
    cout << " numberOfEntries: number of keys, we recommend at least 100000000" << endl;
    cout << " algorithmId: -1 for all (default), or 0..n to only run this algorithm" << endl;
    cout << " algorithmId: can also be a comma-separated list of non-negative integers" << endl;
    for(auto i : names) {
        cout << "           "<< i.first << " : " << i.second << endl;
    }
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
  int algorithmId = -1;
  std::set<int> algos;
  if (argc > 2) {
      if(strstr(argv[2],",") != NULL) {
        // we have a list of algos
        algorithmId = 9999999; // disabling
        parse_comma_separated(argv[2], algos);
      } else {
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
  vector<uint64_t> to_add = seed == -1 ?
      GenerateRandom64Fast(add_count, rand()) :
      GenerateRandom64Fast(add_count, seed);
  vector<uint64_t> to_lookup = seed == -1 ?
      GenerateRandom64Fast(SAMPLE_SIZE, rand()) :
      GenerateRandom64Fast(SAMPLE_SIZE, seed + add_count);

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

  assert(to_lookup.size() == SAMPLE_SIZE);
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
  size_t actual_sample_size = SAMPLE_SIZE;
  if (actual_sample_size > add_count) {
    cout << "WARNING: Your set contains only " << add_count << ". We can't very well support a sample size of " <<   SAMPLE_SIZE << endl;
    actual_sample_size = add_count;
  }

  if (actual_sample_size > to_lookup.size()) {
    throw out_of_range("to_lookup must contain at least SAMPLE_SIZE values");
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
  if (algorithmId == 0 || algorithmId < 0 || (algos.find(0) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter<uint64_t, uint8_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[0] << cf << endl;
  }

  if (algorithmId == 1 || algorithmId < 0 || (algos.find(1) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter2<uint64_t, uint32_t, UInt12Array, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[1] << cf << endl;
  }

  if (algorithmId == 2 || algorithmId < 0 || (algos.find(2) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter<uint64_t, uint16_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[2] << cf << endl;
  }

  if (algorithmId == 3 || algorithmId < 0 || (algos.find(3) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 8, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[3]<< cf << endl;
  }

  if (algorithmId == 4 || algorithmId < 0 || (algos.find(4) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 12, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[4] << cf << endl;
  }

  if (algorithmId == 5 || algorithmId < 0 || (algos.find(5) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 16, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[5] << cf << endl;
  }

  if (algorithmId == 6 || algorithmId < 0 || (algos.find(6) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 13, PackedTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[6] << cf << endl;
  }

  if (algorithmId == 7 || algorithmId < 0 || (algos.find(7) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 8, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[7] << cf << endl;
  }

  if (algorithmId == 8 || algorithmId < 0 || (algos.find(8) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 12, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[8]<< cf << endl;
  }

  if (algorithmId == 9 || algorithmId < 0 || (algos.find(9) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 16, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[9] << cf << endl;
  }

#ifdef __AVX2__
  if (algorithmId == 10 || algorithmId < 0 || (algos.find(10) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed<SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[10] << cf << endl;
  }
#endif

  if (algorithmId == 11 || algorithmId < 0 || (algos.find(11) != algos.end())) {
      auto start_time = NowNanos();
      std::sort(to_add.begin(), to_add.end());
      const auto sort_time = NowNanos() - start_time;
      std::cout << "Sort time: " << sort_time / to_add.size() << " ns/key\n";
  }

  if (algorithmId == 12 || algorithmId < 0 || (algos.find(12) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilterPlus<uint64_t, uint8_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[12] << cf << endl;
  }

  if (algorithmId == 13 || algorithmId < 0 || (algos.find(13) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilterPlus<uint64_t, uint16_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[13] << cf << endl;
  }

  if (algorithmId == 14 || algorithmId < 0 || (algos.find(14) != algos.end())) {
      auto cf = FilterBenchmark<
          GcsFilter<uint64_t, 8, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[14] << cf << endl;
  }

#ifdef __AVX2__
  if (algorithmId == 15 || algorithmId < 0 || (algos.find(15) != algos.end())) {
      auto cf = FilterBenchmark<
          GQFilter<uint64_t, 8, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << names[15] << cf << endl;
  }
#endif

// other algorithms, but not all that interesting or
// not fully optimized

#ifdef __AVX2__
  if (algorithmId == 16 || (algos.find(16) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilter<SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << "BlockedBloom-2^n" << cf << endl;
  }
#endif

  if (algorithmId == 17 || (algos.find(17) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 8, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << "Cuckoo8-2^n" << cf << endl;
  }

  if (algorithmId == 18 || (algos.find(18) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 12, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << "Cuckoo12-2^n" << cf << endl;
  }

  if (algorithmId == 19 || (algos.find(19) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 16, SingleTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << "Cuckoo16-2^n" << cf << endl;
  }

  if (algorithmId == 20 || (algos.find(20) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 13, PackedTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << "CuckooSemiSort13-2^n" << cf << endl;
  }

  if (algorithmId == 21 || (algos.find(21) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter2n<uint64_t, uint8_t, UIntArray<uint8_t>, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << "Xor8-2^n" << cf << endl;
  }

  if (algorithmId == 22 || (algos.find(22) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter2<uint64_t, uint16_t, NBitArray<uint16_t, 10>, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << "Xor10" << cf << endl;
  }

  if (algorithmId == 23 || (algos.find(23) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter2<uint64_t, uint16_t, NBitArray<uint16_t, 14>, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << "Xor14" << cf << endl;
  }

  if (algorithmId == 24 || (algos.find(24) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter2<uint64_t, uint32_t, UInt10Array, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << "Xor10.x" << cf << endl;
  }

  if (algorithmId == 25 || (algos.find(25) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter10<uint64_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << "Xor10" << cf << endl;
  }

  if (algorithmId == 26 || (algos.find(26) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter10_666<uint64_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << "Xor10.666" << cf << endl;
  }

  if (algorithmId == 27 || (algos.find(27) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter13<uint64_t, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << "Xor13" << cf << endl;
  }

  if (algorithmId == 37 || algorithmId < 0 || (algos.find(37) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 8, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[37] << cf << endl;
  }

  if (algorithmId == 38 || algorithmId < 0 || (algos.find(38) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 12, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[38] << cf << endl;
  }

  if (algorithmId == 39 || algorithmId < 0 || (algos.find(39) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 16, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[39] << cf << endl;
  }

#ifdef __AVX2__
  if (algorithmId == 40 || algorithmId < 0 || (algos.find(40) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed<SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed, true);
      cout << setw(NAME_WIDTH) << names[40] << cf << endl;
  }
#endif



// broken algorithms (don't always find all key)
/*
  if (algorithmId == 25) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 9, PackedTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << "CuckooSemiSort9-2^n" << cf << endl;
  }

  if (algorithmId == 26) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 17, PackedTable, SimpleMixSplit>>(
          add_count, to_add, distinct_add, to_lookup, distinct_lookup, intersectionsize, hasduplicates, mixed_sets, seed);
      cout << setw(NAME_WIDTH) << "CuckooSemiSort17-2^n" << cf << endl;
  }
*/

}
