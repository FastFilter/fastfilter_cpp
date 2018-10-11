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

#include <stdio.h>

#include "cuckoofilter.h"
#include "cuckoofilter_stable.h"
#include "xorfilter.h"
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
  double adds_per_nano;
  map<int, double> finds_per_nano; // The key is the percent of queries that were expected
                                   // to be positive
  double false_positive_probabilty;
  double bits_per_item;
};

// Output for the first row of the table of results. type_width is the maximum number of
// characters of the description of any table type, and find_percent_count is the number
// of different lookup statistics gathered for each table. This function assumes the
// lookup expected positive probabiilties are evenly distributed, with the first being 0%
// and the last 100%.
string StatisticsTableHeader(int type_width, int find_percent_count) {
  ostringstream os;

  os << string(type_width, ' ');
  os << setw(12) << right << "million";
  for (int i = 0; i < find_percent_count; ++i) {
    os << setw(8) << "find";
  }
  os << setw(8) << "" << setw(11) << "" << setw(11)
     << "optimal" << setw(8) << "wasted" << setw(8) << "million" << endl;

  os << string(type_width, ' ');
  os << setw(12) << right << "adds/sec";
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
  constexpr double NANOS_PER_MILLION = 1000;
  os << fixed << setprecision(2) << setw(12) << right
     << stats.adds_per_nano * NANOS_PER_MILLION;
  for (const auto& fps : stats.finds_per_nano) {
    os << setw(8) << fps.second * NANOS_PER_MILLION;
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
  }
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
  }
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
  }
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
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<XorFilter<ItemType, FingerprintType, HashFamily>> {
  using Table = XorFilter<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename FingerprintStorageType, typename HashFamily>
struct FilterAPI<XorFilter2<ItemType, FingerprintType, FingerprintStorageType, HashFamily>> {
  using Table = XorFilter2<ItemType, FingerprintType, FingerprintStorageType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename FingerprintStorageType, typename HashFamily>
struct FilterAPI<XorFilter2n<ItemType, FingerprintType, FingerprintStorageType, HashFamily>> {
  using Table = XorFilter2n<ItemType, FingerprintType, FingerprintStorageType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<XorFilterPlus<ItemType, FingerprintType, HashFamily>> {
  using Table = XorFilterPlus<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<GcsFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = GcsFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
  }
  static void AddAll(const vector<ItemType> keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
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
  }
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
  }
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


template <typename Table>
Statistics FilterBenchmark(
    size_t add_count, const vector<uint64_t>& to_add, const vector<uint64_t>& to_lookup, int seed) {
  if (add_count > to_add.size()) {
    throw out_of_range("to_add must contain at least add_count values");
  }
  size_t actual_sample_size = SAMPLE_SIZE;
  if (actual_sample_size > add_count) {
    cout << "WARNING: Your set contains only " << add_count << ". We can't very well support a sample size of " <<   SAMPLE_SIZE << endl;
    actual_sample_size = add_count;
  }

  if (actual_sample_size > to_lookup.size()) {
    throw out_of_range("to_lookup must contain at least SAMPLE_SIZE values");
  }
  size_t distinct_lookup;
  size_t distinct_add;
  size_t intersectionsize = match_size(to_lookup, to_add, &distinct_lookup, & distinct_add);
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
  Table filter = FilterAPI<Table>::ConstructFromAddCount(add_count);
  Statistics result;

  // Add values until failure or until we run out of values to add:
  auto start_time = NowNanos();

  for (size_t added = 0; added < add_count; ++added) {
    FilterAPI<Table>::Add(to_add[added], &filter);
  }
  // for the XorFilter
  FilterAPI<Table>::AddAll(to_add, 0, add_count, &filter);
  // sanity check:
  for (size_t added = 0; added < add_count; ++added) {
    assert(FilterAPI<Table>::Contain(to_add[added], &filter) == 1);
  }
  result.add_count = add_count;
  result.adds_per_nano = add_count / static_cast<double>(NowNanos() - start_time);
  result.bits_per_item = static_cast<double>(CHAR_BIT * filter.SizeInBytes()) / add_count;
  ::std::random_device random;
  size_t found_count = 0;
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
#endif

  for (const double found_probability : {0.0, 0.25, 0.50, 0.75, 1.00}) {
    uint64_t mixingseed = seed == -1 ? random() : seed;
    const auto to_lookup_mixed = DuplicateFreeMixIn(&to_lookup[0], &to_lookup[actual_sample_size], &to_add[0],
    &to_add[add_count], found_probability, mixingseed);

    if(! hasduplicates ) assert(! has_duplicates(to_lookup_mixed));
    assert(to_lookup_mixed.size() == actual_sample_size);
    size_t true_match = match_size(to_lookup_mixed,to_add, NULL, NULL);
    double trueproba =  true_match /  static_cast<double>(actual_sample_size) ;
    double bestpossiblematch = fabs(round(found_probability * actual_sample_size) / static_cast<double>(actual_sample_size) - found_probability);
    double tolerance = bestpossiblematch > 0.01 ? bestpossiblematch : 0.01;
    double probadiff = fabs(trueproba - found_probability);
    if(probadiff >= tolerance) {
      cerr << "WARNING: You claim to have a find proba. of " << found_probability << " but actual is " << trueproba << endl;
    }
    size_t found_before = found_count;
    const auto start_time = NowNanos();
#ifdef __linux__
    unified.start();
#endif
    for (const auto v : to_lookup_mixed) {
      found_count += FilterAPI<Table>::Contain(v, &filter);
    }
#ifdef __linux__
    unified.end(results);
    printf("cycles = %10.zu (cycles per key %10.3f) instructions = %10.zu (ins/key %10.3f,ins/cycles %10.3f) cache misses = %10.zu (misses per keys %10.3f) branch misses = %10.zu (misses per keys %10.3f) \n",
      (size_t)results[0], results[0]*1.0/to_lookup_mixed.size(), (size_t)results[1], results[1]*1.0/to_lookup_mixed.size() , results[1]*1.0/results[0], (size_t)results[2], results[2]*1.0/to_lookup_mixed.size(),
      (size_t)results[3], results[3] * 1.0/to_lookup_mixed.size());
#endif

    const auto lookup_time = NowNanos() - start_time;
    size_t found_this_section = found_count - found_before;
    if (found_this_section < true_match) {
           cerr << "ERROR: Expected to find at least " << true_match << " found " << found_this_section << endl;
           cerr << "ERROR: This is a potential bug!" << endl;
    }
    if (found_probability == 1.00) {
        if (found_this_section != to_lookup_mixed.size()) {
           cerr << "ERROR: Expected to find " << to_lookup_mixed.size() << " found " << found_this_section << endl;
           cerr << "ERROR: Actual intersection is " << true_match << endl;
           cerr << "ERROR: This is a potential bug!" << endl;
        }
    }
    result.finds_per_nano[100 * found_probability] =
        actual_sample_size / static_cast<double>(lookup_time);
    if (0.0 == found_probability) {
      ////////////////////////////
      // This is obviously technically wrong!!! The assumption is that there is no overlap between the random
      // queries and the random content. This is likely true if your 64-bit values were generated randomly,
      // but not true in general.
      ///////////////////////////
      // result.false_positive_probabilty =
      //    found_count / static_cast<double>(to_lookup_mixed.size());
      if(to_lookup_mixed.size() == intersectionsize) {
        cerr << "WARNING: fpp is probably meaningless! " << endl;
      }
      result.false_positive_probabilty = (found_count  - intersectionsize) / static_cast<double>(to_lookup_mixed.size() - intersectionsize);
    }
  }
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

int main(int argc, char * argv[]) {

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <numberOfEntries> [<algorithmId> [<seed>]]" << endl;
    cerr << " numberOfEntries: number of keys" << endl;
    cerr << " algorithmId: -1 for all (default), or 0..n to only run this algorithm" << endl;
    cerr << " seed: seed for the PRNG; -1 for random seed (default)" << endl;
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
  if (argc > 2) {
      stringstream input_string_2(argv[2]);
      input_string_2 >> algorithmId;
      if (input_string_2.fail()) {
          cerr << "Invalid number: " << argv[2];
          return 2;
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
      GenerateRandom64<::std::random_device>(add_count) :
      GenerateRandom64Fast(add_count, seed);
  vector<uint64_t> to_lookup = seed == -1 ?
      GenerateRandom64<::std::random_device>(SAMPLE_SIZE) :
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

  constexpr int NAME_WIDTH = 32;

  cout << StatisticsTableHeader(NAME_WIDTH, 5) << endl;

  if (algorithmId == 0 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          XorFilter<uint64_t, uint8_t, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Xor8" << cf << endl;
  }

  if (algorithmId == 1 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          XorFilter2<uint64_t, uint32_t, UInt12Array, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Xor12" << cf << endl;
  }

  if (algorithmId == 2 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          XorFilter<uint64_t, uint16_t, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Xor16" << cf << endl;
  }

  if (algorithmId == 3 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 8, SingleTable, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "CuckooStable8" << cf << endl;
  }

  if (algorithmId == 4 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 12, SingleTable, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "CuckooStable12" << cf << endl;
  }

  if (algorithmId == 5 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 16, SingleTable, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "CuckooStable16" << cf << endl;
  }

  if (algorithmId == 6 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 13, PackedTable, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "CuckooSemiSortStable13" << cf << endl;
  }

  if (algorithmId == 7 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 8, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Bloom8" << cf << endl;
  }

  if (algorithmId == 8 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 12, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Bloom12" << cf << endl;
  }

  if (algorithmId == 9 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 16, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Bloom16" << cf << endl;
  }

#ifdef __AVX2__
  if (algorithmId == 10 || algorithmId < 0) {
      auto cf = FilterBenchmark<SimdBlockFilter<SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "BlockedBloom" << cf << endl;
  }
#endif

  if (algorithmId == 11 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          GcsFilter<uint64_t, 8, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "GCS" << cf << endl;
  }

#ifdef __AVX2__
  if (algorithmId == 12 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          GQFilter<uint64_t, 8, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "CQF" << cf << endl;
  }
#endif

  if (algorithmId == 13 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          XorFilterPlus<uint64_t, uint8_t, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Xor+8" << cf << endl;
  }

  if (algorithmId == 14 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          XorFilterPlus<uint64_t, uint16_t, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Xor+16" << cf << endl;
  }

#ifdef __AVX2__
  if (algorithmId == 15 || algorithmId < 0) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed<SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "BlockedBloomFixed" << cf << endl;
  }
#endif

  if (algorithmId == 16 || algorithmId < 0) {
      auto start_time = NowNanos();
      std::sort(to_add.begin(), to_add.end());
      const auto sort_time = NowNanos() - start_time;
      std::cout << "Sort time: " << sort_time / to_add.size() << " ns/key\n";
  }

// other algorithms, but not all that interesting or
// not fully optimized
  if (algorithmId == 17 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 8, SingleTable, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Cuckoo2^n-8" << cf << endl;
  }

  if (algorithmId == 18 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 12, SingleTable, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Cuckoo2^n-12" << cf << endl;
  }

  if (algorithmId == 19 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 16, SingleTable, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Cuckoo2^n-16" << cf << endl;
  }

  if (algorithmId == 20 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 13, PackedTable, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "CuckooSemiSort2^n-13" << cf << endl;
  }

  if (algorithmId == 21 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          XorFilter2n<uint64_t, uint8_t, UIntArray<uint8_t>, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Xor2^n-8" << cf << endl;
  }

  if (algorithmId == 22 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          XorFilter2<uint64_t, uint16_t, NBitArray<uint16_t, 10>, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Xor-10" << cf << endl;
  }

  if (algorithmId == 23 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          XorFilter2<uint64_t, uint16_t, NBitArray<uint16_t, 14>, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "Xor-14" << cf << endl;
  }


// broken algorithms (don't always find all key)
/*
  if (algorithmId == 24 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 9, PackedTable, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "SemiSort-9-2^n" << cf << endl;
  }

  if (algorithmId == 25 || algorithmId < 0) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 17, PackedTable, SimpleMixSplit>>(
          add_count, to_add, to_lookup, seed);
      cout << setw(NAME_WIDTH) << "SemiSort-17-2^n" << cf << endl;
  }
*/

}
