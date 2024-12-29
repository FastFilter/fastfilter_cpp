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
#include "random.h"
#include "timing.h"
#ifdef __linux__
#include "linux-perf-events.h"
#endif
#include "filterapi.h"

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


// Output for the first row of the table of results. type_width is the maximum number of
// characters of the description of any table type, and find_percent_count is the number
// of different lookup statistics gathered for each table. This function assumes the
// lookup expected positive probabiilties are evenly distributed, with the first being 0%
// and the last 100%.
string StatisticsTableHeader(int type_width, const std::vector<double> &found_probabilities) {
  ostringstream os;

  os << string(type_width, ' ');
  os << setw(8) << right << "";
  os << setw(8) << right << "";
  for (size_t i = 0; i < found_probabilities.size(); ++i) {
    os << setw(8) << "find";
  }
  os << setw(8) << "1*add+";
  os << setw(8) << "" << setw(11) << "" << setw(11)
     << "optimal" << setw(8) << "wasted" << setw(8) << "million" << endl;

  os << string(type_width, ' ');
  os << setw(8) << right << "add";
  os << setw(8) << right << "remove";
  for (double prob : found_probabilities) {
    os << setw(8 - 1) << static_cast<int>(prob * 100.0) << '%';
  }
  os << setw(8) << "3*find";
  os << setw(9) << "ε%" << setw(11) << "bits/item" << setw(11)
     << "bits/item" << setw(8) << "space%" << setw(8) << "keys";
  return os.str();
}

// Overloading the usual operator<< as used in "std::cout << foo", but for Statistics
template <class CharT, class Traits>
basic_ostream<CharT, Traits>& operator<<(
    basic_ostream<CharT, Traits>& os, const Statistics& stats) {
  os << fixed << setprecision(2) << setw(8) << right
     << stats.nanos_per_add;
  double add_and_find = 0;
  os << fixed << setprecision(2) << setw(8) << right
     << stats.nanos_per_remove;
  for (const auto& fps : stats.nanos_per_finds) {
    os << setw(8) << fps.second;
    add_and_find += fps.second;
  }
  add_and_find = add_and_find * 3 / stats.nanos_per_finds.size();
  add_and_find += stats.nanos_per_add;
  os << setw(8) << add_and_find;

  // we get some nonsensical result for very small fpps
  if(stats.false_positive_probabilty > 0.0000001) {
    const auto minbits = log2(1 / stats.false_positive_probabilty);
    os << setw(8) << setprecision(4) << stats.false_positive_probabilty * 100
       << setw(11) << setprecision(2) << stats.bits_per_item << setw(11) << minbits
       << setw(8) << setprecision(1) << 100 * (stats.bits_per_item / minbits - 1)
       << " " << setw(7) << setprecision(3) << (stats.add_count / 1000000.);
  } else {
    os << setw(8) << setprecision(4) << stats.false_positive_probabilty * 100
       << setw(11) << setprecision(2) << stats.bits_per_item << setw(11) << 64
       << setw(8) << setprecision(1) << 0
       << " " << setw(7) << setprecision(3) << (stats.add_count / 1000000.);
  }
  return os;
}


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
    size_t add_count, vector<uint64_t> to_add,
    size_t intersectionsize, const std::vector<samples_t> & mixed_sets, bool batchedadd = false, bool remove = false) {
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
  printf("cycles: %5.1f/key, instructions: (%5.1f/key, %4.2f/cycle) cache misses: %5.2f/key branch misses: %4.2f/key effective frequency %4.2f GHz\n",
    results[0]*1.0/add_count,
    results[1]*1.0/add_count ,
    results[1]*1.0/results[0],
    results[2]*1.0/add_count,
    results[3]*1.0/add_count,
    results[0]*1.0/time);
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
    printf("cycles: %5.1f/key, instructions: (%5.1f/key, %4.2f/cycle) cache misses: %5.2f/key branch misses: %4.2f/key effective frequency %4.2f GHz\n",
      results[0]*1.0/to_lookup_mixed.size(),
      results[1]*1.0/to_lookup_mixed.size(),
      results[1]*1.0/results[0],
      results[2]*1.0/to_lookup_mixed.size(),
      results[3] * 1.0/to_lookup_mixed.size(),
      results[0]*1.0/lookup_time);
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
    printf("cycles: %5.1f/key, instructions: (%5.1f/key, %4.2f/cycle) cache misses: %5.2f/key branch misses: %4.2f/key effective frequency %4.2f GHz\n",
      results[0]*1.0/add_count,
      results[1]*1.0/add_count ,
      results[1]*1.0/results[0],
      results[2]*1.0/add_count,
      results[3]*1.0/add_count,
      results[0]*1.0/time);
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
        if (ss.peek() == ',') {
            ss.ignore();
        } else if (ss.peek() == '-') {
            ss.ignore();
            int j;
            ss >> j;
            for (i++; i <= j; i++) {
                answer.insert(i);
            }
        }
    }
}


int main(int argc, char * argv[]) {
  std::map<int,std::string> names = {
    // Xor
    {1000, "Xor8-naive"}, {1002, "Xor16-naive"},
    {0, "Xor8"}, {2, "Xor16"},
    {3, "Xor+8"}, {4, "Xor+16"},
    // Cuckooo
    {10,"Cuckoo8"}, {11,"Cuckoo12"}, {12,"Cuckoo16"},
    {13,"CuckooSemiSort13"},
    {14, "Cuckoo8-2^n"}, {15, "Cuckoo12-2^n"}, {16, "Cuckoo16-2^n"},
    {17, "CuckooSemiSort13-2^n"},
    {18, "CuckooFuse8"},
    {19, "CuckooFuse16"},
    // GCS
    {20,"GCS"},
#ifdef __AVX2__
    // CQF + VQF
    {30,"CQF"},
    {31,"VQF"},
    // TwoChoicer
    {32,"TwoChoice"},
    // Prefix
    {35,"PF[TC]"},
    {36,"PF[CF-12-Flex]"},
    {37,"PF[BBF-Flex]"},
#endif
    // Bloom
    {40, "Bloom8"}, {41, "Bloom12" }, {42, "Bloom16"},
    {43, "Bloom8-addAll"}, {44, "Bloom12-addAll"}, {45, "Bloom16-addAll"},
    {46, "BranchlessBloom8-addAll"},
    {47, "BranchlessBloom12-addAll"},
    {48, "BranchlessBloom16-addAll"},
    // Blocked Bloom
    {50, "SimpleBlockedBloom"},
#ifdef __aarch64__
    {51, "BlockedBloom"},
    {52, "BlockedBloom-addAll"},
#elif defined( __AVX2__)
    {51, "BlockedBloom"},
    {52, "BlockedBloom-addAll"},
    {53, "BlockedBloom64"},
#endif
#ifdef __SSE4_1__
    {54, "BlockedBloom16"},
#endif

    // Counting Bloom
    {60, "CountingBloom10-addAll"},
    {61, "SuccCountingBloom10-addAll"},
    {62, "SuccCountBlockBloom10"},
    {63, "SuccCountBlockBloomRank10"},

    {70, "Xor8-singleheader"},
    {72, "BinaryFuse8-singleheader"},

    {80, "Morton"},

    {96, "XorBinaryFuse8-naive"},
    {97, "XorBinaryFuse16-naive"},
    {116, "XorBinaryFuse8"},
    {117, "XorBinaryFuse16"},
    {118, "XorBinaryFuse8-4wise"},
    {119, "XorBinaryFuse16-4wise"},
    {1056, "HomogRibbon64_5"},
    {1076, "HomogRibbon64_7"}, // interesting
    {1086, "HomogRibbon64_8"},
    {1096, "HomogRibbon64_9"},
    {1115, "HomogRibbon32_11"},
    {1116, "HomogRibbon64_11"},
    {1136, "HomogRibbon64_13"},
    {1156, "HomogRibbon64_15"},
    {1776, "HomogRibbon64_7.7"},
    {2056, "BalancedRibbon64Pack_5"},
    {2076, "BalancedRibbon64Pack_7"},
    {2086, "BalancedRibbon64Pack_8"},
    {2096, "BalancedRibbon64Pack_9"},
    {2116, "BalancedRibbon64Pack_11"},
    {2136, "BalancedRibbon64Pack_13"},
    {2156, "BalancedRibbon64Pack_15"},
    {2776, "BalancedRibbon64Pack_7.7"},
    {3056, "StandardRibbon64_5"},
    {3072, "StandardRibbon64_25PctPad_7"},
    {3073, "StandardRibbon64_20PctPad_7"},
    {3074, "StandardRibbon64_15PctPad_7"},
    {3075, "StandardRibbon64_10PctPad_7"},
    {3076, "StandardRibbon64_7"},
    {3086, "StandardRibbon64_8"},
    {3088, "StandardRibbon64_8_Smash"},
    {3096, "StandardRibbon64_9"},
    {3116, "StandardRibbon64_11"},
    {3136, "StandardRibbon64_13"},
    {3156, "StandardRibbon64_15"},
    {3776, "StandardRibbon64_7.7"},

    // Sort
    {100, "Sort"},
  };

  // Parameter Parsing ----------------------------------------------------------
  const char * add_count_str;

  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <numberOfEntries> [<algorithmId> [<seed>]]" << endl;
    cout << " numberOfEntries: number of keys, we recommend at least 100000000" << endl;
    cout << " algorithmId: -1 for all default algos, or 0..n to only run this algorithm" << endl;
    cout << " algorithmId: can also be a comma-separated list of non-negative integers, or ranges (e.g. 1,10-19,30)" << endl;
    for(auto i : names) {
      cout << "           "<< i.first << " : " << i.second << endl;
    }
    cout << " algorithmId: can also be set to the string 'all' if you want to run them all, including some that are excluded by default" << endl;
    cout << " seed: seed for the PRNG; -1 for random seed (default)" << endl;
    cout << endl;
    add_count_str = "10000000";  
  } else {
    add_count_str = argv[1];
  }
  stringstream input_string(add_count_str);
  size_t add_count;
  input_string >> add_count;
  if (input_string.fail()) {
    cerr << "Invalid number: " << add_count_str << endl;
    return 2;
  }
  int algorithmId = -1; // -1 is just the default
  std::set<int> algos;
  if (argc > 2) {
      if(strcmp(argv[2],"all") == 0) {
         for(auto i : names) {// we add all the named algos.
           algos.insert(i.first);
         }
      } else if(strstr(argv[2],",") != NULL || 
          (strstr(argv[2],"-") != NULL && argv[2][0] != '-')) {
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

  if(intersectionsize > 0) {
    cout << "WARNING: Out of the lookup table, "<< intersectionsize<< " ("<<intersectionsize * 100.0 / to_lookup.size() << "%) of values are present in the filter." << endl;
  }

  if(distinct_lookup != to_lookup.size()) {
    cout << "WARNING: Lookup contains "<< (to_lookup.size() - distinct_lookup)<<" duplicates." << endl;
  }
  if(distinct_add != to_add.size()) {
    cout << "WARNING: Filter contains "<< (to_add.size() - distinct_add) << " duplicates." << endl;
  }

  if (actual_sample_size > to_lookup.size()) {
    std::cerr << "actual_sample_size = "<< actual_sample_size << std::endl;
    throw out_of_range("to_lookup must contain at least actual_sample_size values");
  }

  std::vector<samples_t> mixed_sets;

  const std::vector<double> found_probabilities = { 0.0, 0.25, 0.5, 0.75, 1.0 };

  for (const double found_probability : found_probabilities) {
    std::cout << "generating samples with probability " << found_probability <<" ... " << std::flush;

    struct samples thisone;
    thisone.found_probability = found_probability;
    thisone.actual_sample_size = actual_sample_size;
    uint64_t mixingseed = seed == -1 ? random() : seed;
    // seed could be 0 (incremental numbers, or random() might return 0), which we can't use
    if (seed == 0) seed = 1;
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
  cout << StatisticsTableHeader(NAME_WIDTH, found_probabilities) << endl;

  // Algorithms ----------------------------------------------------------
  int a;

  // Xor ----------------------------------------------------------
  a = 0;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter<uint64_t, uint8_t, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 2;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilter<uint64_t, uint16_t, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 1000;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          naive::XorFilter<uint64_t, uint8_t, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  a = 1002;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          naive::XorFilter<uint64_t, uint16_t, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 2000;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          prefetch::XorFilter<uint64_t, uint8_t, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  a = 2002;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          prefetch::XorFilter<uint64_t, uint16_t, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilterPlus<uint64_t, uint8_t, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 4;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorFilterPlus<uint64_t, uint16_t, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  // Cuckoo ----------------------------------------------------------
  a = 10;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 8, SingleTable, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 11;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 12, SingleTable, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 12;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 16, SingleTable, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 13;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilterStable<uint64_t, 13, PackedTable, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 14;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 8, SingleTable, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 15;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 12, SingleTable, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 16;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 16, SingleTable, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 17;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFilter<uint64_t, 13, PackedTable, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 18;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFuseFilter<uint64_t, uint8_t>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 19;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CuckooFuseFilter<uint64_t, uint16_t>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  // GCS ----------------------------------------------------------
  a = 20;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          GcsFilter<uint64_t, 8, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  // CQF ----------------------------------------------------------
#ifdef __AVX2__
  a = 30;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          GQFilter<uint64_t, 8, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 31;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
     auto cf = FilterBenchmark<
          VQFilter<uint64_t, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true, false);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
#endif
#if __PF_AVX512__
  a = 32;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
            TC_shortcut<SimpleMixSplit>>(
            add_count, to_add, intersectionsize, mixed_sets, false, false /* set to true to support deletions. */); 
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  // Prefix ---------------------------------------------------------
    a = 35;
    if (algorithmId == a || (algos.find(a) != algos.end())) {
        auto cf = FilterBenchmark<
                Prefix_Filter<TC_shortcut<SimpleMixSplit>>>(
                add_count, to_add, intersectionsize, mixed_sets,  false, false);
        cout << setw(NAME_WIDTH) << names[a] << cf << endl;
    }
    a = 36;
    if (algorithmId == a || (algos.find(a) != algos.end())) {
        auto cf = FilterBenchmark<
                Prefix_Filter<CuckooFilterStable<uint64_t, 12, SingleTable, SimpleMixSplit>>>(
                add_count, to_add, intersectionsize, mixed_sets,  false, false);
        cout << setw(NAME_WIDTH) << names[a] << cf << endl;
    }
    a = 37;
    if (algorithmId == a || (algos.find(a) != algos.end())) {
        auto cf = FilterBenchmark<
                Prefix_Filter<SimdBlockFilterFixed<SimpleMixSplit>>>(
                add_count, to_add, intersectionsize, mixed_sets,  false, false);
        cout << setw(NAME_WIDTH) << names[a] << cf << endl;
    }
#endif

  // Bloom ----------------------------------------------------------
  a = 40;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 8, false, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 41;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 12, false, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 42;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 16, false, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 43;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 8, false, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 44;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 12, false, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 45;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 16, false, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 46;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 8, true, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 47;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 12, true, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 48;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 16, true, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  a = 48;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BloomFilter<uint64_t, 16, true, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  // Blocked Bloom ----------------------------------------------------------
  a = 50;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          SimpleBlockFilter<8, 8, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
#ifdef __aarch64__
  a = 51;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed<SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 52;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed<SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
#endif
#ifdef __AVX2__
  a = 51;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed<SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 52;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed<SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 53;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
        auto cf = FilterBenchmark<SimdBlockFilterFixed64<SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
#endif
#ifdef __SSE4_1__
  a = 54;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<SimdBlockFilterFixed16<SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
#endif

  // Counting Bloom ----------------------------------------------------------
  a = 60;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          CountingBloomFilter<uint64_t, 10, true, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 61;
  if (algorithmId == a  || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          SuccinctCountingBloomFilter<uint64_t, 10, true, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  true, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 62;
  if (algorithmId == a  || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          SuccinctCountingBlockedBloomFilter<uint64_t, 10, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 63;
  if (algorithmId == a  || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          SuccinctCountingBlockedBloomRankFilter<uint64_t, 10, SimpleMixSplit>>(
          add_count, to_add, intersectionsize, mixed_sets,  false, true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  a = 70;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          XorSingle>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  a = 72;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BinaryFuseSingle>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  a = 80;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          MortonFilter>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  // Xor Binary Fuse Filter ----------------------------------------------------------
  a = 96;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          xorbinaryfusefilter_naive::XorBinaryFuseFilter<uint64_t, uint8_t>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 97;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          xorbinaryfusefilter_naive::XorBinaryFuseFilter<uint64_t, uint16_t>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 102;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          xorbinaryfusefilter_naive4wise::XorBinaryFuseFilter<uint64_t, uint8_t>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 103;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          xorbinaryfusefilter_naive4wise::XorBinaryFuseFilter<uint64_t, uint16_t>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  // Xor Binary Fuse Filter ----------------------------------------------------------
  a = 116;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          xorbinaryfusefilter_lowmem::XorBinaryFuseFilter<uint64_t, uint8_t>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 117;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          xorbinaryfusefilter_lowmem::XorBinaryFuseFilter<uint64_t, uint16_t>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 118;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          xorbinaryfusefilter_lowmem4wise::XorBinaryFuseFilter<uint64_t, uint8_t>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 119;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          xorbinaryfusefilter_lowmem4wise::XorBinaryFuseFilter<uint64_t, uint16_t>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }  
  
   // Homogeneous Ribbon
  a = 1056;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          HomogRibbonFilter<uint64_t, 5>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 1076;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          HomogRibbonFilter<uint64_t, 7>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 1086;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          HomogRibbonFilter<uint64_t, 8>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 1096;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          HomogRibbonFilter<uint64_t, 9>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 1116;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          HomogRibbonFilter<uint64_t, 11>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 1136;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          HomogRibbonFilter<uint64_t, 13>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 1156;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          HomogRibbonFilter<uint64_t, 15>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  a = 2056;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BalancedRibbonFilter<uint64_t, 5, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 2076;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BalancedRibbonFilter<uint64_t, 7, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 2086;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BalancedRibbonFilter<uint64_t, 8, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 2096;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BalancedRibbonFilter<uint64_t, 9, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 2116;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BalancedRibbonFilter<uint64_t, 11, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 2136;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BalancedRibbonFilter<uint64_t, 13, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 2156;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BalancedRibbonFilter<uint64_t, 15, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 2776;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          BalancedRibbonFilter<uint64_t, 0, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }

  a = 3056;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 5, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3072;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 7, 25>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3073;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 7, 20>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3074;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 7, 15>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3075;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 7, 10>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3076;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 7, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3086;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 7, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3088;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 7, 0, true>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3096;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 9, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3116;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 11, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3136;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 13, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3156;
  if (algorithmId == a || algorithmId < 0 || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 15, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
      cout << setw(NAME_WIDTH) << names[a] << cf << endl;
  }
  a = 3776;
  if (algorithmId == a || (algos.find(a) != algos.end())) {
      auto cf = FilterBenchmark<
          StandardRibbonFilter<uint64_t, 0, 0>>(
          add_count, to_add, intersectionsize, mixed_sets,  true);
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
  if(to_add.size() < 100000) {
      std::cout << "You specified a relatively small input size; we recommend running the benchmark multiple times in such cases.\n";
  }


}
