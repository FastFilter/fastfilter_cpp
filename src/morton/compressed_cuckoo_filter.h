/*
Copyright (c) 2019 Advanced Micro Devices, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Author: Alex D. Breslow
        Advanced Micro Devices, Inc.
        AMD Research
*/
#ifndef _COMPRESSED_CUCKOO_FILTER_H
#define _COMPRESSED_CUCKOO_FILTER_H

#include <bitset>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <type_traits> // For std::conditional
#include <set>
#include <unordered_set>
#include <algorithm>

#include "fixed_point.h"
#include "block.h"
#include "morton_util.h"
#include "hash_util.h"
#include "vector_types.h"
#include "compressed_cuckoo_config.h"
#include "bf.h"

#ifndef INLINE
#define INLINE __attribute__((always_inline)) inline
#endif

struct Tester; // Forward declaration
std::ostream& operator<<(std::ostream& os, __uint128_t integer);

namespace CompressedCuckoo{
  struct AccessCounter{
      uint8_t primary_count = 0;
      uint8_t block_overflow_count = 0;
      uint8_t bucket_overflow_count = 0;
      uint8_t block_and_bucket_overflow_count = 0;
  };

  // **** Morton filter ****
  // t_ota_len_bits is the length of the overflow tracking array in bits.
  // If the value is 0, then the functionality should default to a compressed
  // cuckoo filter.
  // I use heavy templating to ensure compile-time definition of class
  // parameters.  I wanted to have a constexpr constructor, but since not
  // every member of the class can be compile-time constant, I was stuck with
  // good old templates.  If you think there's a better solution, let me know.
  // There are some examples of how to instantiate these template parameters in
  // morton_sample_configs.h.
  template<
    uint16_t t_slots_per_bucket,
    uint16_t t_fingerprint_len_bits,
    uint16_t t_ota_len_bits,
    uint16_t t_block_size_bits,
    SerializedFixedPoint t_target_compression_ratio,
    CounterReadMethodEnum t_read_counters_method,
    FingerprintReadMethodEnum t_read_fingerprints_method,
    ReductionMethodEnum t_reduction_method,
    AlternateBucketSelectionMethodEnum t_alternate_bucket_selection_method,
    OverflowTrackingArrayHashingMethodEnum t_morton_ota_hashing_method,
    bool t_resizing_enabled,
    bool t_remap_enabled,
    bool t_collision_resolution_enabled,
    bool t_morton_filter_functionality_enabled,
    bool t_block_fullness_array_enabled,
    bool t_handle_conflicts,
    FingerprintComparisonMethodEnum t_fingerprint_comparison_method
  >
  struct CompressedCuckooFilter{
    //FIXME: Recomment back in later private: // Default but explicitly state
    using block_t =  Morton::Block<t_block_size_bits, t_fingerprint_len_bits,
      atom_t>;
    constexpr static bool _DEBUG = false; // For debug prints
    // Set to true when running measure_bucket_accesses
    constexpr static bool _print_access_counts = false;
    constexpr static bool _resizing_enabled = t_resizing_enabled;
    constexpr static bool _use_bloom_ota = false; // The bit vector was faster.

    // Reserve the zero fingerprint for marking empties.  Although we don't
    // explicitly store empty fingerprints from the logical interpretation,
    // this feature is useful for checking the
    // last slot of FSA.  If it is the zero fingerprint, then we know that
    // the block isn't completely full and this saves us an additive reduction.
    // At low to moderately high loads, setting this flag to true
    // improves insertion throughput by about 50 percent.  However, it has a high price.
    // Currently, it doesn't play nice with my biased insertion algorithm at high loads
    // and slightly lowers lookup throughput.  I've thus set it to false since
    // most workloads are read-heavy, but you might want to set it to true
    // depending on your use case.
    constexpr static bool _special_null_fingerprint = false;

    constexpr static uint_fast16_t _slots_per_bucket = t_slots_per_bucket;
    constexpr static uint_fast16_t _fingerprint_len_bits = t_fingerprint_len_bits;
    // Overflow tracking array's length in bits
    constexpr static uint_fast16_t _ota_len_bits = t_ota_len_bits;
    constexpr static uint_fast16_t _block_size_bits = t_block_size_bits;

    // How to perform insertions (first-fit (FIRST_FIT) -- place fingerprint
    // in first bucket with a spare slot, two-choice (TWO_CHOICE) -- check the
    // loads both both candidate blocks and place the fingerprint in the bucket
    // in the block with less load, hybrid (HYBRID_PIECEWISE) -- start off by
    // using first-fit but transition to using two-choice once you reach a
    // threshold).  FIRST_FIT_OPT accesses a single bucket/block up to a
    // threshold, but once that threshold is reached, it fetches both
    // blocks since the probability of a trivial insertion in the first
    // block decreases as the LF increases, enough so to eventually warrant just
    // grabbing both and calling it a day.  HYBRID_SIMPLE inserts using
    // first fit for k-1 batches and then switches to two-choice for one batch,
    // and repeats.  The parameter k is tunable. You probably don't want to
    // touch this.
    constexpr static InsertionMethodEnum _insertion_method =
      InsertionMethodEnum::FIRST_FIT;//HYBRID_SIMPLE;//FIRST_FIT;FIRST_FIT_OPT;//TWO_CHOICE;//HYBRID_PIECEWISE;

    // How fullness counters are read (w/ or w/o supporting reading counters
    // that cross atoms (machine words) in the
    // block store (READ_RAW, READ_SIMPLE, READ_CROSS))
    // For example, if the fullness counter array has a bit offset of zero
    // within the block, and each fullness counter is
    // two bits in length, then READ_SIMPLE will work.  This is because
    // each counter will be fully contained within a single atom.  However, if
    // you use 3-bit counters, READ_SIMPLE can only be used if the counters
    // fall within a single atom (e.g., bits 0 to 62 inclusive when atom_t
    // is a 64-bit unsigned integer and we have 21 3-bit fullness counters).
    // READ_CROSS will trigger reading a single atom twice if the
    // counter appears in a single atom but will read two atoms if
    // the fullness counter is split across two atoms.
    // READ_RAW is an optimization for when fullness counters are always in
    // the zeroth atom of the block and exclusively that atom.
    constexpr static CounterReadMethodEnum _read_counters_method =
      t_read_counters_method;
    // How fingerprints are read (w/ or w/o supporting reading fingerprints
    // that cross atom boundaries (READ_SIMPLE and READ_CROSS supported))
    // READ_BYTE is a special implementation for 8-bit fingerprints that
    // are properly aligned to byte boundaries, but it doesn't appear
    // that much faster than READ_SIMPLE.
    constexpr static FingerprintReadMethodEnum _read_fingerprints_method =
      t_read_fingerprints_method;
    // Options are NAIVE_FULL_EXCLUSIVE_SCAN and POP_CNT
    constexpr static ReductionMethodEnum _reduction_method =
      t_reduction_method;
    // Valid values are TABLE_BASED_OFFSET and FUNCTION_BASED_OFFSET
    constexpr static AlternateBucketSelectionMethodEnum
      _alternate_bucket_selection_method = t_alternate_bucket_selection_method;

    // Specific to Morton filter
    // Performance tip: set to RAW_BUCKET_HASH when _buckets_per_block is a power
    // of 2.  CLUSTERED_BUCKET_HASH is recommended if you have a lot of negative
    // lookups. LEMIRE_FINGERPRINT_MULTIPLY is recommended if you have a number
    // of buckets per block that is smaller than the length of the OTA in bits.
    // At present, I am mostly using CLUSTERED_BUCKET_HASH because it plays
    // nicely with biasing fingerprint kickouts from a subset of buckets.
    constexpr static OverflowTrackingArrayHashingMethodEnum
      _morton_ota_hashing_method = t_morton_ota_hashing_method;

    // Disable remapping (false) if you want raw speed and don't care if some
    // elements can't be stored in the filter.  Enable with true.
    constexpr static bool _remap_enabled = t_remap_enabled;
    constexpr static bool _collision_resolution_enabled =
      t_collision_resolution_enabled;
    constexpr static bool _morton_filter_functionality_enabled =
      (t_ota_len_bits > 0);

    // Set to true to enable the block fullness array, which for each block
    // allocates a single bit that is 0 if the block has any free
    // capacity and 1 if it doesn't (all slots full).
    constexpr static bool _block_fullness_array_enabled =
      t_block_fullness_array_enabled;

    // Detect and manage conflicts on insertions
    constexpr static bool _handle_conflicts = t_handle_conflicts;

    // This enum controls how fingerprints are read from a bucket.
    // VARIABLE_COUNT uses the fullness counter as a guide to skip over
    // logical slots that are empty.
    // FIXED_COUNT_AGGRESSIVE loops over all the logical slots in the
    // the bucket but masks those that are invalid.  It can end up reading
    // past the end of the block, so you have to allocate additional padding
    // to mitigate segmentation faults.
    // SEMI_FIXED is similar to FIXED_COUNT_AGGRESSIVE, but it strives to
    // prevent the reading past the end of filter's storage and thus
    // doesn't require extra memory.
    constexpr static FingerprintComparisonMethodEnum
      _fingerprint_comparison_method = t_fingerprint_comparison_method;

    constexpr static uint_fast8_t _max_pop_count_width_in_bits = 128;

    // Add back in the const if you don't need resizability.
    uint_fast64_t _total_buckets;
    uint_fast64_t _total_slots;

    constexpr static uint_fast64_t _fullness_counter_width =
      util::log2ceil(_slots_per_bucket + 1);
    constexpr static uint_fast64_t _bits_per_uncompressed_bucket =
      t_fingerprint_len_bits * _slots_per_bucket;
    constexpr static uint_fast64_t _buckets_per_block = (_block_size_bits -
      _ota_len_bits)/ (_fullness_counter_width +
      FixedPoint(t_target_compression_ratio).to_double() *
      _bits_per_uncompressed_bucket);
;
    // Used for controlling which overflows require updating the OTA
    // This number is the highest block-local bucket index for which fingerprints
    // that overflow are not recorded in the OTA.  Setting this value above -1
    // means that some overflows will not be recorded and thus on negative
    // lookups for those logical buckets, both candidate buckets will need to
    // be checked.  This is mainly a performance optimization for heavy loads
    // to reduce cruft build up in the OTA that causes performance to
    // drop.  Here, you know right away by examining the block-local bucket
    // index whether you'll likely only need to access one bucket or perhaps
    // two.  Such information can inform prefetching and which lookup or
    // insertion algorithm to use.
    constexpr static int_fast16_t _ota_lbi_insertion_threshold = -1;


    // Add back in const if you don't need resizability
    uint_fast64_t _total_blocks;
    constexpr static uint_fast64_t _max_fingerprints_per_block = (_block_size_bits - _ota_len_bits - _buckets_per_block * _fullness_counter_width) / _fingerprint_len_bits;
    constexpr static uint_fast64_t _fullness_counters_offset = 0;
    constexpr static uint_fast64_t _overflow_tracking_array_offset = _fullness_counters_offset + _buckets_per_block * _fullness_counter_width;
    //const uint_fast64_t _block_full_tracker_offset;
    constexpr static uint_fast64_t _fingerprint_offset = _fullness_counters_offset + _buckets_per_block * _fullness_counter_width + _ota_len_bits;

    atom_t _popcount_masks[max_fullness_counter_width] = {};
    __uint128_t _popcount_masks128[max_fullness_counter_width] = {};

    using fca_t = typename std::conditional<(_fullness_counter_width * _buckets_per_block <= 64), uint64_t, __uint128_t>::type;
    fca_t _reduction_masks[util::log2ceil(_buckets_per_block)] = {};

    counter_t* _summed_counters; //[_buckets_per_block + 1];
    block_t* _storage;
    std::vector<bool> _block_fullness_array;
    BitMixMurmur _hasher; // Yields more consistent performance
    // The number of times that we've doubled the filter's capacity
    uint_fast16_t _resize_count;

    size_t sizeInBytes;

    friend Tester; // Class with a bunch of test functions in test.cc

  public:
  // Constructor
  explicit CompressedCuckooFilter(uint64_t total_slots) :
    // Hashing mechanism requires even number of total buckets
    // We round up to a number of buckets that's even but also for which
    // _total_buckets % _buckets_per_block is 0. So in the worst case,
    // we allocate a little less than two extra blocks worth of storage.
    _total_buckets(determine_total_buckets<uint_fast64_t>(_slots_per_bucket,
      total_slots / FixedPoint(t_target_compression_ratio).to_double(), _buckets_per_block)),
    _total_slots(_total_buckets * _slots_per_bucket), // Logical slots not physical
    _total_blocks(_total_buckets / _buckets_per_block),
    _block_fullness_array(_block_fullness_array_enabled ? _total_blocks : 0, 0),
    _resize_count(0)
  {

    sizeInBytes = 0;
    sizeInBytes += sizeof(bool) * (_block_fullness_array_enabled ? _total_blocks : 0);
    sizeInBytes += sizeof(counter_t) * (_buckets_per_block + 1);
    sizeInBytes += sizeof(block_t) * _total_blocks;

    // Supporting dual use as a compressed cuckoo filter and Morton filter
    static_assert((_morton_filter_functionality_enabled ^ (_ota_len_bits == 0)),
    "ERROR: If Morton filter functionality is enabled, then the overflow tracking array must be at least one bit in length.");

    /* TODO: Get this working.  These checks are important.
    if(_fullness_counter_width * _buckets_per_block > 8 * sizeof(atom_t)){
      // You may only make RAW reads if the fullness counter array is less than
      // or equal to the width of a block store atom.
      static_assert((_read_counters_method != CounterReadMethodEnum::READ_RAW),
        "ERROR: You can only do raw reads if the fullness counter array is "
        " within the first atom of the block store.");

      // You may only make SIMPLE reads if the counter width is a power of 2.
      // FIXME: Include checks for fullness counter array offset.  For instance,
      // _fullness_counters_offset makes it so that the counters no longer
      // evenly divide into an atom, that will also break the code.
      if(_read_counters_method == CounterReadMethodEnum::READ_SIMPLE){
        // Is it a power of 2?
        static_assert((__builtin_popcountll(_fullness_counter_width) == 1),
          "ERROR: Fullness counter width is not a power of 2 when the FCA is "
          " bigger than one atom in size.");
      }
    }
    */

    if(_total_buckets & 0x1){
      std::cerr << "ERROR: Total buckets needs to be even\n";
      exit(1);
    }

    // Check if we can use a popcount-based reduction
    // TODO: If it's possible to know these values at compile time, then
    // perform this check with static_assert.
    if(_reduction_method == ReductionMethodEnum::POP_CNT && !(_fullness_counter_width * _buckets_per_block <= _max_pop_count_width_in_bits)){
      std::cerr << "ERROR: With popcount-based reduction, the fullness counters"
        << " must be less than or equal to the width of a single atom.\n";
      std::cerr << *this << std::endl;
      exit(1);
    }

    // Generate masks for doing exclusive reductions
    generate_popcount_masks<atom_t>(_popcount_masks);
    generate_popcount_masks<__uint128_t>(_popcount_masks128);
    generate_reduction_masks<fca_t>(_reduction_masks);

    // Allocate heap memory so that it's cache aligned.
    heap_allocate_table_and_summed_counters_buffer();
  }

  ~CompressedCuckooFilter(){
    if(g_cache_aligned_allocate){
      free(_summed_counters);
      free(_storage);
    }
    else{
      delete[] _summed_counters;
      delete[] _storage;
    }
  }

  size_t SizeInBytes() {
    return sizeInBytes;
  }

  INLINE atom_t fingerprint_function(const hash_t raw_hash) const{
    atom_t fingerprint =
      raw_hash >> (8 * sizeof(raw_hash) - _fingerprint_len_bits);
    if(_special_null_fingerprint && (fingerprint == 0)){
      fingerprint = 1; // See comment above regarding _special_null_fingerprint
    }
    return fingerprint;
  }

  // Takes the portion of the hash that isn't part of the fingerprint and
  // computes an operation that maps it to a bucket within the table
  INLINE hash_t map_to_bucket(const hash_t raw_hash, const hash_t modulus) const{
    constexpr hash_t one = 1;
    constexpr hash_t shift = 8 * sizeof(raw_hash) - _fingerprint_len_bits;
    constexpr hash_t mask = (one << shift) - one;
    if (_resizing_enabled){
      const hash_t original_modulus = modulus >> _resize_count;
      const hash_t original_bucket_id = util::fast_mod_alternative<hash_t>(
        raw_hash & mask, original_modulus, shift);
      const hash_t original_block_id = original_bucket_id / _buckets_per_block;
      const hash_t lbi = original_bucket_id % _buckets_per_block;
      const atom_t fingerprint = fingerprint_function(raw_hash);
      const hash_t new_block_id = (original_block_id << _resize_count) |
        ((fingerprint >> (_fingerprint_len_bits - _resize_count)) & ((one << _resize_count) - one));
      const hash_t new_bucket_id = (new_block_id * _buckets_per_block) + lbi;
      return new_bucket_id;
    }
    return util::fast_mod_alternative<hash_t>(raw_hash & mask, modulus, shift);
  }

  // The number of buckets in the table must be a power of 2 to use this.
  INLINE hash_t fan_et_al_partial_key_cuckoo_hash_alternate_bucket(hash_t
    bucket_id, const atom_t fingerprint) const{
    return (bucket_id ^ raw_primary_hash(fingerprint)) & (_total_buckets - 1);
  }

  // See comment above
  INLINE hash_t determine_alternate_bucket(hash_t bucket_id,
    const atom_t fingerprint) const{
    int64_t offset;
    constexpr hash_t one = 1;
    switch(_alternate_bucket_selection_method){
      case AlternateBucketSelectionMethodEnum::TABLE_BASED_OFFSET:{
        constexpr int16_t offsets[] = {83, 149, 211, 277, 337, 397, 457, 521,
          587, 653, 719, 787, 853, 919, 983, 1051, 1117, 1181, 1249, 1319, 1399,
          1459,
          1511, 1571, 1637, 1699, 1759, 1823, 1889, 1951, 2017, 1579};//, 1579
        static_assert(offsets[0] > _buckets_per_block,
          "Cannot use TABLE_BASED_OFFSET with so many buckets per block");
        offset = offsets[fingerprint % (sizeof(offsets) / sizeof(offsets[0]))];
        break;
      }

      case AlternateBucketSelectionMethodEnum::FUNCTION_BASED_OFFSET:{
        offset = ((raw_primary_hash(fingerprint) & 0x1fff) +
          (_buckets_per_block)) | one;
        break;
      }
      case AlternateBucketSelectionMethodEnum::FAN_ET_AL_PARTIAL_KEY:
        return fan_et_al_partial_key_cuckoo_hash_alternate_bucket(bucket_id,
          fingerprint);
        break;
    }

    if(_resizing_enabled){
      hash_t block_id = bucket_id / _buckets_per_block;
      hash_t lbi = bucket_id % _buckets_per_block;
      hash_t old_block_id = block_id >> _resize_count;
      hash_t old_bucket_id = old_block_id * _buckets_per_block + lbi;
      bucket_id = old_bucket_id;
    }

    // lbi is synonymous with counter_index in other parts of the code
    offset = (bucket_id & one) ? offset : -offset;

    /* An alternate implementation that works on big tables */
    int64_t output = static_cast<int64_t>(bucket_id) + offset;
    if(_resizing_enabled){
      if(output < 0){
        output += _total_buckets >> _resize_count;
      }
      if(static_cast<uint64_t>(output) >= (_total_buckets >> _resize_count)){
        output -= _total_buckets >> _resize_count;
      }
      // Scale the output to account for resizing
      hash_t output_block_id = output / _buckets_per_block;
      hash_t output_lbi = output % _buckets_per_block;
      hash_t new_output_block_id = (output_block_id << _resize_count) |
        ((fingerprint >> (_fingerprint_len_bits - _resize_count)) & ((one << _resize_count) - one));
      output = new_output_block_id * _buckets_per_block + output_lbi;
    }
    else{
      if(output < 0){
        output += _total_buckets;
      }
      if(static_cast<uint64_t>(output) >= _total_buckets){
        output -= _total_buckets;
      }
    }
    return output;

  }

  inline block_t* allocate_cache_aligned_storage(uint64_t total_blocks){
    size_t allocation_size = sizeof(block_t) * total_blocks;

    block_t* storage;
    const int malloc_failed =
        posix_memalign(reinterpret_cast<void**>(&storage), g_cache_line_size_bytes, allocation_size);
    if (malloc_failed) throw ::std::bad_alloc();

    //block_t* storage = static_cast<block_t*>(aligned_alloc(
    //  g_cache_line_size_bytes, allocation_size));



    // Currently set to false because clear_swath hasn't been rigorously tested
    constexpr bool _only_clear_ota_and_fca = false;
    if(!_only_clear_ota_and_fca){ // Competitive with the code in the loop below

// from https://phabricator.kde.org/D13900
#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && (((__GNUC__ * 100) + __GNUC_MINOR__) >= 800)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif
      memset(storage, 0x0, allocation_size);
#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && (((__GNUC__ * 100) + __GNUC_MINOR__) >= 800)
#pragma GCC diagnostic pop
#endif

      return storage;
    } // ELSE
    for(hash_t i = 0; i < total_blocks; i++){
      // This code assumes that the OTA and FCA appear one after the other
      constexpr uint64_t fca_length_bits = _buckets_per_block *
        _fullness_counter_width;
      constexpr uint64_t combined_fca_ota_len_in_bits = fca_length_bits +
        _ota_len_bits;
      static_assert((_fullness_counters_offset + fca_length_bits ==
        _overflow_tracking_array_offset) && (_fullness_counters_offset <
        _overflow_tracking_array_offset), "In compressed_cuckoo_filter.h,\n"
        "allocate_cache_aligned_storage assumes that the OTA appears directly"
        "\nafter the FCA.  If this is not the case, then you'll have to\n"
        "reimplement this part.");
      constexpr bool use_memset =
        (((fca_length_bits + _ota_len_bits) % 8) == 0)
        && (_fullness_counters_offset % 8 == 0);
      if(use_memset){
        memset(reinterpret_cast<uint8_t*>(&storage[i]) +
          (_fullness_counters_offset / 8), 0x0,
          combined_fca_ota_len_in_bits / 8);
      }
      else{ // TODO: Test me rigorously
        storage[i].clear_swath(_fullness_counters_offset,
          combined_fca_ota_len_in_bits);
      }
    }
    return storage;
  }

  void heap_allocate_table_and_summed_counters_buffer(){
    if(g_cache_aligned_allocate){ // Allocate heap memory so that it's cache
                                // aligned
      // Allocate memory for the block store
      _storage = allocate_cache_aligned_storage(_total_blocks);

      size_t allocation_size = sizeof(*_summed_counters) *
        (_buckets_per_block + 1);

      // _summed_counters = static_cast<decltype(_summed_counters)>(aligned_alloc(g_cache_line_size_bytes, allocation_size));
      const int malloc_failed =
        posix_memalign(reinterpret_cast<void**>(&_summed_counters), g_cache_line_size_bytes, allocation_size);
      if (malloc_failed) throw ::std::bad_alloc();

     // Memset not required, initialized by full_exclusive_scan

      if(_storage == NULL || _summed_counters == NULL){
         std::cerr << "ERROR: Allocating table memory failed" << std::endl;
         exit(1);
      }
    }
    else{
      _summed_counters = new counter_t[_buckets_per_block + 1]();
      _storage = new block_t[_total_blocks]();
    }
  }

  template<class T>
  inline void generate_reduction_masks(T* reduction_masks){
    // TODO: Fix iteration so that it is over the maximum buckets that could
    // appear in a block not the _buckets_per_block
    for(uint32_t i = 0; i < util::log2ceil(_buckets_per_block); i++){
      for(uint32_t j = 0; j < sizeof(T) * 8; j+=2*(_fullness_counter_width << i)){
        reduction_masks[i] +=
          ((static_cast<T>(1) << (_fullness_counter_width << i)) - 1) << j;
      }
      if(_DEBUG){
        std::cout << "REDUCTION MASK " << i << ": " <<
          util::bin_string<T>(reduction_masks[i], _fullness_counter_width)
          << std::endl;
      }
    }
  }

  // Generate masks for using popcount to perform reductions on counters that
  // are small enough to fit multiple per word.  Often, these counters are
  // only several bits wide.  If each counter is 2 bits wide, then we only
  // have to invoke popcount twice, once with the mask 0x55555555 and then with
  // the mask 0xaaaaaaaa.  This code autogenerates those masks for any range of
  // counter widths.
  template <class T>
  inline void generate_popcount_masks(T* popcount_masks){
    for(uint32_t i = 0; i < _fullness_counter_width; i++){
      for(uint32_t j = i; j < sizeof(T) * 8; j+= _fullness_counter_width){
        popcount_masks[i] += (static_cast<T>(1) << j);
      }
      if(_DEBUG){
        std::cout << "POPCOUNT MASK " << i << " " <<
          util::bin_string(popcount_masks[i], _fullness_counter_width)
          << std::endl;
      }
    }
  }

  // It would be better to use an algorithm with logarithmic height that
  // takes advantage of SIMD.  This implementation is naive and has a
  // loop-carried dependency on _summed_counters.
  // Current implementation goes one past a traditional exclusive scan by
  // appending an extra cell that holds the reduction of the entire counter
  // array
  inline counter_t* full_exclusive_scan(const uint64_t block_id){
    _summed_counters[0] = 0;
    for(counter_t i = 1; i < _buckets_per_block + 1; i++){
      _summed_counters[i] = _summed_counters[i - 1] + read_counter(block_id,
        i - 1);
    }
    return _summed_counters;
  }

  // Call the correct implementation depending on which is necessary
  INLINE uint16_t exclusive_reduce_with_popcount(const block_t& b,
    uint8_t counter_index) const{
    constexpr bool use_64bit_popcnt = _fullness_counter_width *
      _buckets_per_block <= 64;
    if(use_64bit_popcnt){
      return exclusive_reduce_with_popcount64(b, counter_index);
    }
    else{
      return exclusive_reduce_with_popcount128(b, counter_index);
    }
  }

  INLINE uint16_t exclusive_reduce(const block_t& b,
    uint8_t counter_index) const{
    uint16_t sum;
    switch(_reduction_method){
      case ReductionMethodEnum::NAIVE_FULL_EXCLUSIVE_SCAN:
        sum = reduce_up_to_index(b, counter_index);
        break;
      case ReductionMethodEnum::POP_CNT:
        sum = exclusive_reduce_with_popcount(b, counter_index);
        break;
      case ReductionMethodEnum::PARALLEL_REDUCE:
        sum = exclusive_reduce_with_parallel_sum(b, counter_index);
        break;
      default:
        std::cerr << "ERROR: Unsupported ReductionMethodEnum value\n";
        exit(1);
        break;
    }
    return sum;
  }

  template<uint64_t SIZE>
  INLINE std::array<counter_t, SIZE> exclusive_reduce_many(
    const std::array<hash_t, SIZE>& block_ids,
    const std::array<counter_t, SIZE>& counter_indexes) const{
    std::array<counter_t, SIZE> sums;
    switch(_reduction_method){
      // Not actually a full scan but a naive acummulator loop
      case ReductionMethodEnum::NAIVE_FULL_EXCLUSIVE_SCAN:
        for(uint64_t i = 0; i < SIZE; i++){
          sums[i] = reduce_up_to_index(_storage[block_ids[i]],
            counter_indexes[i]);
        }
        break;
      case ReductionMethodEnum::POP_CNT:{
        // FIXME: Should ideally check whether the start and end are in the
        // same 64-bit word not the 0th 64-bit word of the block
        constexpr bool use_popcnt64 = (_fullness_counters_offset +
          _fullness_counter_width * _buckets_per_block) <= 64;
        if(use_popcnt64){
          return exclusive_reduce_with_popcount64<batch_size>(block_ids,
            counter_indexes);
        }
        else{
          return exclusive_reduce_with_popcount128<batch_size>(block_ids,
            counter_indexes);
        }
        break;
      }
      case ReductionMethodEnum::PARALLEL_REDUCE:
        for(uint64_t i = 0; i < SIZE; i++){
          sums[i] = exclusive_reduce_with_parallel_sum(_storage[block_ids[i]],
            counter_indexes[i]);
        }
        break;
      default:
        std::cerr << "ERROR: Unsupported ReductionMethodEnum value\n";
        exit(1);
        break;
    }
    return sums;
  }

  // A tree-based reduction scheme for fullness counters to determine where a
  // bucket with block-local logical bucket index counter_index's first
  // fingerprint is stored in the FSA.  A series of masks is applied to ensure
  // that accumulated counts do not overflow by zeroing out fields that
  // have already been tabulated as part of the reduction.  This zeroing
  // is necessary since typical fullness counter widths are one to four bits,
  // which is not enough to count the dozens of fingerprints that are typical
  // for a block that is sized to a hardware cache line.
  uint16_t exclusive_reduce_with_parallel_sum(const block_t& b,
    uint8_t counter_index) const{
    // Mask grabs all counters up to but excluding the counter
    // at counter_index.  Add 1 to counter_index in the
    // product to get the inclusive reduction.
    const fca_t mask = (static_cast<fca_t>(one) << (_fullness_counter_width * counter_index)) - one;
    // FIXME: Assumes the counters are all in the block's least significant bytes
    // Replace with an appropriate read_many_cross as necessary
    // Add necessary safety checks
    fca_t sum;
    memcpy(&sum, &b, sizeof(fca_t));
    sum >>= _fullness_counters_offset;
    sum &= mask;

    // Only apply a mask during the iterations where it is necessary to prevent
    // overflows.
    constexpr uint_fast8_t num_passes = util::log2ceil(_buckets_per_block);
    static_assert(_max_fingerprints_per_block > 0, "The number of fingerprints "
      "per block must be one or more.");
    constexpr uint_fast8_t masked_count = static_cast<uint_fast8_t>(util::log2ceil(
      _max_fingerprints_per_block) / _fullness_counter_width) - 1;
    //constexpr uint_fast8_t masked_count = static_cast<uint_fast8_t>(ceil(log2(
    //  _max_fingerprints_per_block) / _fullness_counter_width)) - 1;
    for(int8_t i = 0; i < masked_count; i++){ // Masked to avoid overflows
      sum = (sum & _reduction_masks[i]) + ((sum >> (_fullness_counter_width << i)) & _reduction_masks[i]);
    }
    for(int8_t i = masked_count; i < num_passes; i++){ // No mask, no overflows
      sum += (sum >> (_fullness_counter_width << i));
    }
    return sum & ((one << (_fullness_counter_width * (masked_count + one))) - one);
  }

  template<uint64_t SIZE>
  INLINE std::array<counter_t, SIZE> exclusive_reduce_with_popcount64(
    const std::array<hash_t, SIZE>& block_ids,
    const std::array<counter_t, SIZE>& counter_indexes) const{
    constexpr uint64_t one = 1;
    uint64_t counters;
    std::array<counter_t, SIZE> sums;

    uint64_t popcount_mask = _popcount_masks[0];
    for(uint64_t i = 0; i < SIZE; i++){
      memcpy(&counters, &_storage[block_ids[i]], sizeof(__uint64_t));
      uint64_t mask = (one << (_fullness_counter_width * (counter_indexes[i])))
        - one;
      counters &= mask;
      sums[i] = 0;
      for(uint8_t j = 0; j < _fullness_counter_width; j++){
        sums[i] += __builtin_popcountll(counters & (popcount_mask << j)) << j;
      }
    }
    return sums;
  }

  template<uint64_t SIZE>
  INLINE std::array<counter_t, SIZE> exclusive_reduce_with_popcount128(
    const std::array<hash_t, SIZE>& block_ids,
    const std::array<counter_t, SIZE>& counter_indexes) const{
    constexpr __uint128_t one = 1;
    constexpr uint64_t uint64_max = ~(static_cast<uint64_t>(0));

    std::array<__uint128_t, SIZE> counters;
    std::array<counter_t, SIZE> sums;

    for(uint64_t i = 0; i < SIZE; i++){
      memcpy(&counters[i], &_storage[block_ids[i]], sizeof(__uint128_t));
      __uint128_t mask = (one << (_fullness_counter_width *
        (counter_indexes[i]))) - one;
      counters[i] &= mask;
      sums[i] = 0;
      __uint128_t popcount_mask128 = _popcount_masks128[0];
      for(uint8_t j = 0; j < _fullness_counter_width; j++){
        const uint64_t popcount_mask1 = popcount_mask128 & uint64_max;
        const uint64_t popcount_mask2 = popcount_mask128 >> 64;
        sums[i] += __builtin_popcountll(counters[i] & popcount_mask1) << j;
        sums[i] += __builtin_popcountll((counters[i] >> 64) & popcount_mask2) << j;
        popcount_mask128 <<= 1;
      } // This version that shifts the mask rather than precomputing it
        // at different offsets is a good deal faster.  Register pressure?
    }
    return sums;
  }

  INLINE uint16_t exclusive_reduce_with_popcount128(const block_t& b,
    uint8_t counter_index) const{
    constexpr __uint128_t one = 1;
    const __uint128_t mask = (one << (_fullness_counter_width * counter_index))
      - one;
    uint8_t sum = 0u;
    __uint128_t counters;
    memcpy(&counters, &b, sizeof(__uint128_t));
    counters &= mask;

    constexpr uint64_t uint64_max = ~(static_cast<uint64_t>(0));

    __uint128_t popcount_mask128 = _popcount_masks128[0]; // We shift this one.

    for(uint8_t i = 0; i < _fullness_counter_width; i++){
      // Use __builtin_popcountll if GCC and not on an X86 system?
      const uint64_t popcount_mask1 = popcount_mask128 & uint64_max;
      const uint64_t popcount_mask2 = popcount_mask128 >> 64;
      sum += __builtin_popcountll(counters & popcount_mask1) << i;
      sum += __builtin_popcountll((counters >> 64) & popcount_mask2) << i;
      popcount_mask128 <<= 1;
    } // This version that shifts the mask rather than precomputing it
      // at different offsets is a good deal faster.  Register pressure?
    return sum;
  }

  // Exclusive reduction on counters that are one or more bits
  // wide but which fit into a 64-bit word
  INLINE uint16_t exclusive_reduce_with_popcount64(const block_t& b,
    uint8_t counter_index) const{
    constexpr atom_t one = 1;
    // Mask grabs all counters up to but excluding the counter
    // at counter_index.  Add 1 to counter_index in the
    // product to get the inclusive reduction.
    const atom_t mask = (one << (_fullness_counter_width * counter_index)) - one;
    uint8_t sum = 0u;
    // FIXME: Assumes the counters are all in the first atom
    // Replace with an appropriate read_many_cross as necessary
    // Add necessary safety checks
    atom_t counters = b[0] >> _fullness_counters_offset;
    // Avoid the part of the reduction that we don't care
    // about.
    counters &= mask;
    atom_t popcount_mask = _popcount_masks[0];
    for(uint8_t i = 0; i < _fullness_counter_width; i++){
      sum += __builtin_popcountll(counters & popcount_mask) << i;
      popcount_mask <<= 1;
    } // This version that shifts the mask rather than precomputing it
      // at different offsets is a good deal faster.  Register pressure?
    return sum;
  }

  // This implementation is just functional.  It's not going to lead to high
  // performance.  It's a serial exclusive reduction.
  counter_t reduce_up_to_index(const block_t& b, uint8_t counter_index) const{
    uint8_t index = 0;  // Index that we start reading the bucket at
    for(uint8_t i = 0; i < counter_index; i++){
      index += b.read_cross(_fullness_counters_offset, _fullness_counter_width,
        i);
    }
    return index;
  }

  INLINE void decrement_fullness_counter(block_t& b, uint32_t index,
    atom_t counter){
    set_fullness_counter(b, index, counter - 1);
  }

  // Overwrites the fullness counter at the index
  INLINE void set_fullness_counter(block_t& b, uint32_t index,
    atom_t counter){
    switch(_read_counters_method){
      case CounterReadMethodEnum::READ_RAW128:
        b.template add_t<__uint128_t>(_fullness_counters_offset,
          _fullness_counter_width, index, counter);
        break;
      // Fullness counter array only in first atom
      case CounterReadMethodEnum::READ_RAW:
        b.add(_fullness_counters_offset, _fullness_counter_width, index,
          counter);
        break;
      case CounterReadMethodEnum::READ_SIMPLE:
        b.add(_fullness_counters_offset, _fullness_counter_width, index,
          counter);
        break;
      case CounterReadMethodEnum::READ_CROSS:
        b.add_cross(_fullness_counters_offset, _fullness_counter_width, index,
          counter);
        break;
    }
  }


  // Increments the fullness counter by 1 at index counter_index in block b
  // See inc_cross for whether it checks for counter overflows
  INLINE void increment_fullness_counter(block_t& b, uint32_t index){
    uint32_t read;
    if(_DEBUG){ // TODO: Replace with debug
      read = b.read_cross(_fullness_counters_offset, _fullness_counter_width,
        index);
    }
    switch(_read_counters_method){
      case CounterReadMethodEnum::READ_CROSS:
        b.inc_cross(_fullness_counters_offset, _fullness_counter_width, index);
        break;
      default: // Individual counters are never split across multiple words
        b.inc(_fullness_counters_offset, _fullness_counter_width, index);
        break;
    }
    if(_DEBUG){
      uint32_t read_after_inc = b.read_cross(_fullness_counters_offset, _fullness_counter_width, index);
      if(read + 1 != read_after_inc){
        std::cout << std::flush;
        std::cerr << "Index: " << index << ". I started with a count of "
          << read <<
          " and ended with a count of " << read_after_inc << std::endl;
        std::cerr << b << std::endl;
      }
    }
  }

  inline hash_t raw_primary_hash(keys_t key) const{
    return _hasher(key);
  }

  inline bool insert_many(const std::vector<keys_t>& keys,
    std::vector<bool>& status, const uint64_t num_keys){
    for(hash_t i = 0; i < num_keys; i += batch_size){
      ar_hash bucket_hashes;
      ar_atom fingerprints;
      for(hash_t j = 0; j < batch_size; j++){
        bucket_hashes[j] = raw_primary_hash(keys[i + j]);
      }
      for(hash_t j = 0; j < batch_size; j++){
        // Now primary buckets
        fingerprints[j] = fingerprint_function(bucket_hashes[j]);
        bucket_hashes[j] = map_to_bucket(bucket_hashes[j], _total_buckets);
      }
      switch(_insertion_method){
        case InsertionMethodEnum::TWO_CHOICE:
          table_store_many_two_choice(bucket_hashes, fingerprints, status, i);
          break;
        case InsertionMethodEnum::FIRST_FIT:
          table_store_many(bucket_hashes, fingerprints, status, i);
          break;
        // Differs from hybrid approach in only several lines of code
        case InsertionMethodEnum::HYBRID_PIECEWISE:
        case InsertionMethodEnum::FIRST_FIT_OPT:{
          double cutoff_fraction = 0.7;//0.89;
          hash_t cutoff_point = cutoff_fraction * _total_blocks *
            _max_fingerprints_per_block;
          if(i < cutoff_point){
            table_store_many(bucket_hashes, fingerprints, status, i);
          }
          else{
            table_store_many_two_choice(bucket_hashes, fingerprints, status, i);
          }
          break;
        }
        case InsertionMethodEnum::HYBRID_SIMPLE:{
          // Execute else for every 1 / cutoff_divisor batches
          constexpr hash_t cutoff_divisor = 3;
          if((i % (cutoff_divisor * batch_size)) != 0){
            table_store_many(bucket_hashes, fingerprints, status, i);
          }
          else{
            table_store_many_two_choice(bucket_hashes, fingerprints, status, i);
          }
          break;
        }
        default: // Put here to make the compiler happy
          std::cerr << "SOMETHING IS WRONG IF YOU ARE HERE\n";
          exit(1);
          break;
      }
    }
    return true; // FIXME: Check statuses
  }

  // Item at a time
  inline bool insert(const keys_t key){
    hash_t raw_hash = raw_primary_hash(key);

    atom_t fingerprint = fingerprint_function(raw_hash);
    // Primary bucket
    hash_t primary_bucket = map_to_bucket(raw_hash, _total_buckets);

    if(_DEBUG){  // Replace with debug variable
      std::cout << "Inserting key: "
        << key << " "
        << "Fingerprint: " << fingerprint << " Bucket index: "
        << primary_bucket << std::endl;
      fflush(stdout);
    }

    uint8_t* counters;
    if(_DEBUG){
      std::cout << "IN DEBUG" << std::endl;
      counters = new uint8_t[_total_buckets];
      for(uint64_t i = 0; i < _total_buckets; i++){
        counters[i] = read_counter(i / _buckets_per_block,
          i % _buckets_per_block);
      }
    }
    bool ret = table_store(primary_bucket, fingerprint);
    if(_DEBUG){
      for(uint64_t i = 0; i < _total_buckets; i++){
        atom_t counter_read = read_counter(i / _buckets_per_block, i % _buckets_per_block);
        if(counters[i] != counter_read && i != primary_bucket){
          std::cerr << "table_store is corrupting counters.\n" << std::endl;
          std::cerr << std::flush;
          exit(1);
        }
        else if(i == primary_bucket && counters[i] + 1u != counter_read){
          std::cerr << "table_store is not properly incrementing counters.\n"
            << std::endl;
          std::cerr << std::flush;
          exit(1);
        }
      }
      delete[] counters;
    }
    return ret;
  }

  inline void likely_contains_many(const std::vector<keys_t>& keys,
    std::vector<bool>& status, const uint64_t num_keys) const{
    for(hash_t i = 0; i < num_keys; i += batch_size){
      ar_hash bucket_hashes;
      ar_atom fingerprints;
      for(hash_t j = 0; j < batch_size; j++){
        bucket_hashes[j] = raw_primary_hash(keys[i + j]);
      }
      for(hash_t j = 0; j < batch_size; j++){
        // Now primary buckets
        fingerprints[j] = fingerprint_function(bucket_hashes[j]);
        bucket_hashes[j] = map_to_bucket(bucket_hashes[j],
          _total_buckets);
      }
      // Write the output statuses directly to the output status vector "status"
      table_read_and_compare_many(bucket_hashes, fingerprints, status, i);
    }
  }

  inline bool table_delete_item(hash_t bucket_id, atom_t fingerprint){
    uint64_t block_id = bucket_id / _buckets_per_block;
    uint16_t counter_index = (bucket_id % _buckets_per_block);

    // How many slots are full in the current bucket
    counter_t full_slots = read_counter(block_id, counter_index);

    // Skip a bunch of computation if it's not necessary
    if(full_slots == 0){
      return false;
    }

    counter_t bucket_start_index = get_bucket_start_index(block_id,
      counter_index);

    bool success = false;
    uint8_t j = 0;
    do{ // You always execute at least one iteration of the loop if you get here
      if(fingerprint == read_fingerprint(block_id, bucket_start_index + j)){
        delete_fingerprint_right_displace(_storage[block_id], bucket_start_index
          + j);
        decrement_fullness_counter(_storage[block_id], counter_index,
          full_slots);
        success = true;
        break;
      }
    } while(++j != full_slots);

    // -1 because we have already decremented the fullness counter
    if(_block_fullness_array_enabled && success & (get_bucket_start_index(
      block_id, _buckets_per_block) == _max_fingerprints_per_block - 1)){
      _block_fullness_array[block_id] = 0;
    }

    return success;
  }

  inline void delete_many(const std::vector<keys_t>& keys,
    std::vector<bool>& status, const uint64_t num_keys){
    for(hash_t i = 0; i < num_keys; i += batch_size){
      ar_hash bucket_hashes;
      ar_atom fingerprints;
      for(hash_t j = 0; j < batch_size; j++){
        bucket_hashes[j] = raw_primary_hash(keys[i + j]);
        fingerprints[j] = fingerprint_function(bucket_hashes[j]);
      }
      for(hash_t j = 0; j < batch_size; j++){
        // Now primary buckets
        bucket_hashes[j] = map_to_bucket(bucket_hashes[j],
          _total_buckets);
      }
      table_delete_item_many(bucket_hashes, fingerprints, status, i);
    }
  }

  inline void table_delete_item_many(const ar_hash& bucket_ids, const ar_atom&
    fingerprints, std::vector<bool>& status, const hash_t write_offset){
    ar_hash block_ids;
    ar_counter counter_indexes;
    if(_handle_conflicts){
      // Number of buckets in the Bloom filter
      constexpr uint64_t num_buckets = 64; // Must be power of 2
      BlockedBF::BloomFilter<num_buckets> bf;
      std::bitset<batch_size> conflict_vector;
      for(uint_fast32_t i = 0; i < batch_size; i++){
        block_ids[i] = bucket_ids[i] / _buckets_per_block;
        counter_indexes[i] = bucket_ids[i] % _buckets_per_block;
        conflict_vector[i] = conflict_exists<num_buckets>(bf, block_ids[i]);
      }

      // Fall back to one at a time processing if two updates would be
      // applied to the same block in a batch
      if(conflict_vector.any()){
        for(uint_fast64_t i = 0; i < batch_size; i++){
          status[write_offset + i] = table_delete_item(bucket_ids[i],
            fingerprints[i]);
          if(status[write_offset + i] == false){
            hash_t alternate_bucket_id =
              determine_alternate_bucket(bucket_ids[i], fingerprints[i]);
            status[write_offset + i] = table_delete_item(alternate_bucket_id,
              fingerprints[i]);
          }
        }
        return;
      }
    }
    else{ // Not handling conflicts
      for(uint_fast64_t i = 0; i < batch_size; i++){
        block_ids[i] = bucket_ids[i] / _buckets_per_block;
        counter_indexes[i] = bucket_ids[i] % _buckets_per_block;
      }
    }

    ar_counter full_slots;
    read_counter_many(block_ids, counter_indexes, full_slots);

    ar_counter i_with_secondary_lookup;
    counter_t secondary_count = 0;

    ar_counter bucket_start_indexes;
    for(uint_fast64_t i = 0; i < batch_size; i++){
      bucket_start_indexes[i] = get_bucket_start_index(block_ids[i],
        counter_indexes[i]);
      uint8_t discovery_slot = return_slot_id_on_match<>(block_ids[i],
        bucket_start_indexes[i], full_slots[i], fingerprints[i]);

      if(discovery_slot != _slots_per_bucket){
        delete_fingerprint_right_displace(_storage[block_ids[i]],
          bucket_start_indexes[i] + discovery_slot);
        decrement_fullness_counter(_storage[block_ids[i]], counter_indexes[i],
          full_slots[i]);
        status[write_offset + i] = true;
        if(_block_fullness_array_enabled & (get_bucket_start_index(
          block_ids[i], _buckets_per_block) == _max_fingerprints_per_block - 1)){
          _block_fullness_array[block_ids[i]] = 0;
        }
      }
      else{
        i_with_secondary_lookup[secondary_count] = i;
        secondary_count++;
      }
    }
    // This has to be a separate loop.  Otherwise, an additional round of
    // conflict detection is necessary for secondary items since they may
    // throw off the precomputed counts for primary items.
    for(uint_fast64_t i = 0; i < secondary_count; i++){
      uint_fast64_t true_i = i_with_secondary_lookup[i];
      hash_t secondary_bucket_id = determine_alternate_bucket(
        bucket_ids[true_i], fingerprints[true_i]);
      status[write_offset + true_i] =
        table_delete_item(secondary_bucket_id, fingerprints[true_i]);
    }
  }

  // Item at a time
  inline bool delete_item(const keys_t key){
    hash_t raw_hash = raw_primary_hash(key);
    atom_t fingerprint = fingerprint_function(raw_hash);
    // Primary bucket
    hash_t primary_bucket = map_to_bucket(raw_hash, _total_buckets);
    bool return_status = false;
    if(!_remap_enabled){
      return_status = table_delete_item(primary_bucket, fingerprint);
    }
    else{ // Delete from first bucket.  Only proceed to second bucket if
          // a matching fingerprint is not found in the first bucket.
      return_status = table_delete_item(primary_bucket, fingerprint);
      if(!return_status){
        hash_t secondary_bucket = determine_alternate_bucket(primary_bucket,
          fingerprint);
        return_status = table_delete_item(secondary_bucket, fingerprint);
        // TODO: Try to clear the bit that the overflown item had set
        //attempt_to_clear_ota_bit(primary_bucket, secondary_bucket, fingerprint);
      }
    }
    return return_status;
  }

  // Item at a time
  inline bool likely_contains(const keys_t key){
    hash_t raw_hash = raw_primary_hash(key);
    atom_t fingerprint = fingerprint_function(raw_hash);
    // Primary bucket
    hash_t primary_bucket = map_to_bucket(raw_hash, _total_buckets);

    // Idealized implementation with no remapping necessary
    if(!_remap_enabled){
      return table_read_and_compare(primary_bucket, fingerprint);
    }

    // TODO: Check that this still works
    // Compressed cuckoo filter implementation
    else if(!_morton_filter_functionality_enabled){
      hash_t secondary_bucket = determine_alternate_bucket(primary_bucket,
        fingerprint);
      return table_read_and_compare(primary_bucket, fingerprint) |
        table_read_and_compare(secondary_bucket, fingerprint);
    }

    // Morton filter implementation
    else{
      bool status1 = table_read_and_compare(primary_bucket, fingerprint);
      bool status2 = false;
      if(get_overflow_status(primary_bucket, fingerprint)){ // Possible overflow?
        hash_t secondary_bucket = determine_alternate_bucket(primary_bucket,
          fingerprint);
        status2 = table_read_and_compare(secondary_bucket, fingerprint);
      }
      return status1 | status2;
    }
  }

  inline counter_t get_bucket_start_index(uint64_t block_id, uint16_t
    counter_index) const{
    counter_t bucket_start_index;
    switch(_reduction_method){
      case(ReductionMethodEnum::POP_CNT): // TODO: check _fullness_counter_width * _buckets_per_block <= sizeof(atom_t) * 8llu)
        bucket_start_index = exclusive_reduce_with_popcount(_storage[block_id],
          counter_index);
        break;
      // Naive but not full exclusive scan
      case(ReductionMethodEnum::NAIVE_FULL_EXCLUSIVE_SCAN):
        bucket_start_index = reduce_up_to_index(_storage[block_id],
          counter_index);
        break;
      case(ReductionMethodEnum::PARALLEL_REDUCE):
        bucket_start_index = exclusive_reduce_with_parallel_sum(
          _storage[block_id], counter_index);
        break;
      default:
        std::cerr << "ERROR: Undefined setting for _reduction_method\n";
        exit(1);
        break;
    }
    return bucket_start_index;
  }

  // Reports the C parameter from the VLDB'18 paper.  This is the ratio of physical slots in the FSA
  // to logical slots per block.
  constexpr double report_compression_ratio() const {
    return static_cast<double>(_max_fingerprints_per_block) / (_buckets_per_block * _slots_per_bucket);
  }

  // This is the $\alpha_C$ term in the paper.
  inline double report_block_occupancy(){
    uint64_t full_slots_count = 0;
    for(uint64_t block_id = 0; block_id < _total_blocks; block_id++){
      full_slots_count += get_bucket_start_index(block_id, _buckets_per_block);
    }
    return static_cast<double>(full_slots_count) / (_max_fingerprints_per_block * _total_blocks);
  }

  // Methods specific to Morton filters

  // Reports what fraction of the bits of the Overflow Tracking Array are set
  double report_ota_occupancy(){
    uint64_t set_bit_count = 0;
    if(_morton_filter_functionality_enabled){
      for(uint64_t i = 0; i < _total_blocks; i++){
        for(uint64_t bit_index = 0; bit_index < _ota_len_bits; bit_index++){
          set_bit_count += _storage[i].read_bit(_overflow_tracking_array_offset
            + bit_index);
        }
      }
    }
    return static_cast<double>(set_bit_count) / (_total_blocks * _ota_len_bits);
  }

  void test_fingerprint_in_bucket_many_morton(const ar_hash& bucket_ids,
    const ar_hash& block_ids,
    const ar_counter& bucket_start_indexes, const ar_counter& full_slots,
    const ar_atom& fingerprints, std::vector<bool>& status,
    const hash_t write_offset) const{
    for(uint_fast32_t i = 0; i < batch_size; i++){
      bool found_finger = test_fingerprint_in_bucket<>(block_ids[i],
        bucket_start_indexes[i],
        full_slots[i], fingerprints[i]);

      uint_fast8_t secondary_lookup_necessary = (!found_finger) &
        get_overflow_status(bucket_ids[i], fingerprints[i]);
      if(secondary_lookup_necessary){
         hash_t secondary_bucket_id = determine_alternate_bucket(
          bucket_ids[i],
          fingerprints[i]);
        hash_t block_id = secondary_bucket_id / _buckets_per_block;
        uint16_t counter_index = secondary_bucket_id % _buckets_per_block;
        uint16_t bucket_start_index = get_bucket_start_index(block_id,
          counter_index);
        counter_t secondary_full_slots = read_counter(block_id, counter_index);

        found_finger = test_fingerprint_in_bucket<>(
          block_id,
          bucket_start_index,
          secondary_full_slots, fingerprints[i]);
      }
      status[write_offset + i] = found_finger;
    }
  }

  inline void test_fingerprint_in_bucket_many(const ar_hash& block_ids,
    const ar_counter& bucket_start_indexes, const ar_counter& full_slots,
    const ar_atom& fingerprints, std::vector<bool>& status,
    const hash_t write_offset) const{
    if(_DEBUG){
      util::print_array<ar_hash>("block_ids: ", block_ids);
      util::print_array<ar_counter>("bucket_start_indexes: ",
        bucket_start_indexes);
      util::print_array<ar_counter>("full_slots: ", full_slots);
      util::print_array<ar_atom>("fingerprints: ", fingerprints);
      std::cout << "Write Offset: " << write_offset << std::endl;
    }
    for(uint_fast32_t i = 0; i < batch_size; i++){
      status[write_offset + i] = test_fingerprint_in_bucket<>(block_ids[i],
        bucket_start_indexes[i], full_slots[i], fingerprints[i]);
    }
  }

  template<FingerprintComparisonMethodEnum t_comparison_method =
    _fingerprint_comparison_method>
  INLINE bool test_fingerprint_in_bucket(hash_t block_id,
    counter_t bucket_start_index, uint8_t full_slots,
    atom_t fingerprint) const{
    return return_slot_id_on_match<t_comparison_method>(block_id,
      bucket_start_index, full_slots, fingerprint) != _slots_per_bucket;
  }

  // Returns the index of the fingerprint in the bucket.  If it doesn't find
  // it, then it returns _slots_per_bucket.
  template<FingerprintComparisonMethodEnum t_comparison_method =
    _fingerprint_comparison_method>
  INLINE hash_t return_slot_id_on_match(hash_t block_id,
    counter_t bucket_start_index, uint8_t full_slots, atom_t fingerprint) const{

    hash_t match_index = _slots_per_bucket; // Slot ID of match within a bucket

    // I had these series of if statements as a single switch statement,
    // but GCC wasn't properly evaluating the conditional at compile time
    // despite t_comparison_method being known at compile time.
    // I saw performance degradation of about 30%.  Thus, I added the if
    // statements.

    // The main configuration.  I found this to give the best performance
    // on large filters even though the loop count is variable.
    if (t_comparison_method == FingerprintComparisonMethodEnum::VARIABLE_COUNT) {
      // This loop executes a varying number of times.
      // Check below for a version that doesn't.
      for(uint8_t i = 0; i < full_slots; i++){
        if(fingerprint == read_fingerprint(block_id, bucket_start_index + i)){
          match_index = i;
          //break; You would think this would be a good idea but often it is not
        }
      }
    }

    else if (t_comparison_method == FingerprintComparisonMethodEnum::FIXED_COUNT_AGGRESSIVE) {
      // Requires padding the allocated memory for the block store
      // this algorithmic variant can read past the end of a block
      // by several bytes.  Make sure you pad the memory if you use this.
      // No variable loop count
      for(uint8_t i = 0; (i < _slots_per_bucket); i++){
        bool valid_read = (i < full_slots) & (bucket_start_index + i < _max_fingerprints_per_block);
        bool found = (fingerprint == read_fingerprint(block_id,
          bucket_start_index + i)) & valid_read;
        if(found){
          match_index = i;
        }
      }
    }

    // Same as above but stops reading if you would read past of the end
    // of the block.  The approach above instead performs more masking.
    else if (t_comparison_method == FingerprintComparisonMethodEnum::SEMI_FIXED) {
      for(uint8_t i = 0; (i < _slots_per_bucket) &
        (bucket_start_index + i < _max_fingerprints_per_block); i++){
        bool valid_read = i < full_slots;
        bool found = (fingerprint == read_fingerprint(block_id,
          bucket_start_index + i)) & valid_read;
        if(found){
          match_index = i;
        }
      }
    }
    return match_index;
  }

  inline void table_read_and_compare_many_pessimistic(const ar_hash& bucket_ids,
    const ar_atom& fingerprints, std::vector<bool>& status,
    const hash_t write_offset) const{
    ar_hash block_ids[2];
    ar_hash alternate_bucket_ids;
    ar_counter counter_indexes[2];
    for(uint_fast32_t i = 0; i < batch_size; i++){
      alternate_bucket_ids[i] = determine_alternate_bucket(bucket_ids[i],
        fingerprints[i]);
      block_ids[0][i] = bucket_ids[i] / _buckets_per_block;
      block_ids[1][i] = alternate_bucket_ids[i] / _buckets_per_block;
      counter_indexes[0][i] = bucket_ids[i] % _buckets_per_block;
      counter_indexes[1][i] = alternate_bucket_ids[i] % _buckets_per_block;
    }
    ar_counter full_slots[2];
    read_counter_many(block_ids[0], counter_indexes[0], full_slots[0]);
    read_counter_many(block_ids[1], counter_indexes[1], full_slots[1]);

    ar_counter bucket_start_indexes[2];
    for(uint_fast32_t i = 0; i < batch_size; i++){
      bucket_start_indexes[0][i] = get_bucket_start_index(block_ids[0][i],
        counter_indexes[0][i]);
      bucket_start_indexes[1][i] = get_bucket_start_index(block_ids[1][i],
        counter_indexes[1][i]);
    }

    for(uint_fast32_t i = 0; i < batch_size; i++){
      std::bitset<_slots_per_bucket> found1;
      std::bitset<_slots_per_bucket> found2;
      for(uint_fast8_t j = 0; j < full_slots[0][i]; j++){
        found1[j] = (fingerprints[i] == read_fingerprint(block_ids[0][i],
          bucket_start_indexes[0][i] + j));
      }
      for(uint_fast8_t j = 0; j < full_slots[1][i]; j++){
        found2[j] = (fingerprints[i] == read_fingerprint(block_ids[1][i],
          bucket_start_indexes[1][i] + j));
      }
      status[write_offset + i] = found1.any() | found2.any();
    }
  }

  inline void table_read_and_compare_many(const ar_hash& bucket_ids,
    const ar_atom& fingerprints, std::vector<bool>& status,
    const hash_t write_offset) const{
    ar_hash block_ids;
    ar_counter counter_indexes;
    for(uint_fast32_t i = 0; i < batch_size; i++){
      block_ids[i] = bucket_ids[i] / _buckets_per_block;
      counter_indexes[i] = bucket_ids[i] % _buckets_per_block;
    }
    ar_counter full_slots;
    read_counter_many(block_ids, counter_indexes, full_slots);

    ar_counter bucket_start_indexes;
    for(uint_fast32_t i = 0; i < batch_size; i++){
      bucket_start_indexes[i] = get_bucket_start_index(block_ids[i],
        counter_indexes[i]);
    }

    // Idealized implementation with no remapping necessary
    if(!_remap_enabled){
      test_fingerprint_in_bucket_many(block_ids, bucket_start_indexes,
        full_slots, fingerprints, status, write_offset);
    }

    // Compressed cuckoo filter implementation
    else if(!_morton_filter_functionality_enabled){
      // Check first bucket
      test_fingerprint_in_bucket_many(block_ids, bucket_start_indexes,
        full_slots, fingerprints, status, write_offset);

      ar_hash secondary_bucket_ids;
      ar_hash secondary_block_ids;
      ar_counter secondary_counter_indexes;
      ar_counter secondary_bucket_start_indexes;
      ar_counter secondary_full_slots;

      for(uint_fast32_t i = 0; i < batch_size; i++){
        secondary_bucket_ids[i] = determine_alternate_bucket(bucket_ids[i],
          fingerprints[i]);
        secondary_block_ids[i] = secondary_bucket_ids[i] / _buckets_per_block;
        secondary_counter_indexes[i] =
          secondary_bucket_ids[i] % _buckets_per_block;
        secondary_bucket_start_indexes[i] = get_bucket_start_index(
          secondary_block_ids[i], secondary_counter_indexes[i]);
      }
      read_counter_many(secondary_block_ids, secondary_counter_indexes,
        secondary_full_slots);

      // Check second bucket
      test_fingerprint_in_bucket_many(secondary_block_ids,
        secondary_bucket_start_indexes,
        secondary_full_slots, fingerprints, status, write_offset);
    }

    // Morton Filter
    else{
      test_fingerprint_in_bucket_many_morton(bucket_ids, block_ids,
        bucket_start_indexes, full_slots, fingerprints, status, write_offset);
    }

  }

  // Global slot ID is the global index of the first slot of the bucket we're
  // searching.  Function returns true if fingerprint matches an item stored
  // in the bucket.
  inline bool table_read_and_compare(uint64_t bucket_id,
    atom_t fingerprint) const{
    uint64_t block_id = bucket_id / _buckets_per_block;
    uint16_t counter_index = (bucket_id % _buckets_per_block);
    // How many slots are full in the current bucket
    counter_t full_slots = read_counter(block_id, counter_index);

    // Retrieve item
    counter_t bucket_start_index = get_bucket_start_index(block_id,
      counter_index);

    return test_fingerprint_in_bucket<>(block_id, bucket_start_index,
      full_slots, fingerprint);
  }

  // Set full_slots to the number of slots that are full per bucket
  inline void read_counter_many(const ar_hash& block_ids,
    const ar_counter& counter_indexes, ar_counter& full_slots) const{
    for(uint_fast32_t i = 0; i < batch_size; i++){
      full_slots[i] = read_counter(block_ids[i], counter_indexes[i]);
    }
  }

  // In place writing of the fingerprint without left shifting
  INLINE void write_fingerprint(block_t& block, uint64_t fsa_slot_id,
    atom_t fingerprint){
    // Check if we can use a faster method for setting fingerprints in the FSA
    constexpr bool fsa_using_byte_multiples = ((_fingerprint_len_bits % 8 == 0)
      || (_fingerprint_len_bits == 4)) && (_fingerprint_offset % 8 == 0);
    if(fsa_using_byte_multiples){
      block.add(_fingerprint_offset, _fingerprint_len_bits, fsa_slot_id,
        fingerprint);
    }
    else{
      block.add_cross(_fingerprint_offset, _fingerprint_len_bits, fsa_slot_id,
        fingerprint);
    }
  }

  INLINE void write_fingerprint_left_displace(block_t& block,
    uint64_t fsa_slot_id, atom_t fingerprint){
    // Check if we can use a faster method for setting fingerprints in the FSA
    constexpr bool fsa_using_byte_multiples = (_fingerprint_len_bits % 8 == 0)
      && (_fingerprint_offset % 8 == 0);
    if(fsa_using_byte_multiples){
      constexpr uint64_t fsa_size_bytes = (_max_fingerprints_per_block *
        _fingerprint_len_bits) / 8;
      constexpr uint64_t fingerprint_len_bytes = _fingerprint_len_bits / 8;
      constexpr uint64_t fingerprint_offset_bytes = _fingerprint_offset / 8;
      block.add_left_displace(fingerprint_offset_bytes, fingerprint_len_bytes,
        fsa_slot_id, fingerprint, fsa_size_bytes);
    }
    else{
      block.add_cross_left_displace(_fingerprint_offset, _fingerprint_len_bits,
        fsa_slot_id, fingerprint);
    }
  }

  INLINE void delete_fingerprint_right_displace(block_t& block,
    uint64_t fsa_slot_id){
    // Check if we can use a faster method for deleting fingerprints in the FSA
    constexpr bool fsa_using_byte_multiples = (_fingerprint_len_bits % 8 == 0)
      && (_fingerprint_offset % 8 == 0);
    if(fsa_using_byte_multiples){
      constexpr uint64_t fsa_size_bytes = (_max_fingerprints_per_block *
        _fingerprint_len_bits) / 8;
      constexpr uint64_t fingerprint_len_bytes = _fingerprint_len_bits / 8;
      constexpr uint64_t fingerprint_offset_bytes = _fingerprint_offset / 8;
      block.del_right_displace(fingerprint_offset_bytes, fingerprint_len_bytes,
        fsa_slot_id, fsa_size_bytes);
    }
    else{
      block.del_cross_right_displace(_fingerprint_offset, _fingerprint_len_bits,
        fsa_slot_id);
    }
  }

  // The block_slot_id tells you which fingerprint to select within the block
  // E.g., if there is capacity for 35 fingerprints, block_slot_id can vary
  // from 0 to 34 inclusive.
  INLINE atom_t read_fingerprint(uint64_t block_id, uint64_t block_slot_id)
    const{
    atom_t fingerprint = 0;
    switch(_read_fingerprints_method){
      // Accessing up to 2 atoms of the block store
      case FingerprintReadMethodEnum::READ_CROSS:
        fingerprint = _storage[block_id].read_cross(_fingerprint_offset,
          _fingerprint_len_bits, block_slot_id);
        break;
      case FingerprintReadMethodEnum::READ_SIMPLE:
        fingerprint = _storage[block_id].read(_fingerprint_offset,
          _fingerprint_len_bits, block_slot_id);
        break;
      case FingerprintReadMethodEnum::READ_BYTE:{
        // Get the index of the byte that we want to read
        uint64_t byte_index = (_fingerprint_offset / 8) + block_slot_id;
        fingerprint = _storage[block_id].read_byte(byte_index);
        break;
      }
      default:
        std::cerr << "ERROR: Unknown fingerprint reading mode" << std::endl;
        exit(1);
        break;
    }
    return fingerprint;
  }

  // Method abstracts away which implementation of reading that we use on
  // the block store
  INLINE counter_t read_counter(uint64_t block_id, uint64_t counter_index)
    const{
    counter_t counter = 0;
    switch(_read_counters_method){ // Compile-time constant
      // Accessing up to 2 atoms of the block store
      case CounterReadMethodEnum::READ_CROSS:
        counter = _storage[block_id].read_cross(_fullness_counters_offset,
          _fullness_counter_width, counter_index);
        break;
      // Accessing at most 1 atom of the block store
      case CounterReadMethodEnum::READ_SIMPLE:
        counter = _storage[block_id].read(_fullness_counters_offset,
          _fullness_counter_width, counter_index);
        break;
      // If it's always in atom 0 at offset 0, then just read it
      case CounterReadMethodEnum::READ_RAW:
        counter = _storage[block_id].read_atom0(_fullness_counters_offset,
          _fullness_counter_width, counter_index);
        break;
      case CounterReadMethodEnum::READ_RAW128:
        counter = _storage[block_id].template read_zeroth_word<__uint128_t>(
          _fullness_counters_offset, _fullness_counter_width,
          counter_index);
        break;
      default:
        std::cerr << "ERROR: Unknown counter reading mode" << std::endl;
        exit(1);
        break;
    }
    return counter;
  }

  inline void two_choice_store_many(const ar_hash& bucket_ids_1,
    const ar_hash& bucket_ids_2, const ar_atom& fingerprints,
    const ar_hash& block_ids_1, const ar_hash& block_ids_2,
    std::vector<bool>& statuses, const hash_t offset){
    ar_hash elements_in_blocks_1, elements_in_blocks_2;
    for(uint32_t i = 0; i < batch_size; i++){
      elements_in_blocks_1[i] = exclusive_reduce(
        _storage[block_ids_1[i]], _buckets_per_block);
      elements_in_blocks_2[i] = exclusive_reduce(
        _storage[block_ids_2[i]], _buckets_per_block);
    }
    std::bitset<batch_size> try_first_block_insert;
    for(uint32_t i = 0; i < batch_size; i++){
      switch(_insertion_method){
        case InsertionMethodEnum::FIRST_FIT_OPT:
          try_first_block_insert[i] = (elements_in_blocks_1[i] !=
            _max_fingerprints_per_block) | (elements_in_blocks_2[i] ==
            _max_fingerprints_per_block);
          break;
        default: // HYBRID_PIECEWISE, TWO_CHOICE
          try_first_block_insert[i] = (elements_in_blocks_1[i] <=
            elements_in_blocks_2[i]);
          break;
      }
    }
    // Perfom a blend
    ar_hash bucket_ids;
    ar_hash block_ids;
    ar_counter counter_indexes;
    ar_counter elements_in_blocks;
    for(uint32_t i = 0; i < batch_size; i++){
      bucket_ids[i] = (try_first_block_insert[i] ?
        bucket_ids_1[i] : bucket_ids_2[i]);
      block_ids[i] = (try_first_block_insert[i] ?
        block_ids_1[i] : block_ids_2[i]);
      elements_in_blocks[i] = (try_first_block_insert[i] ?
        elements_in_blocks_1[i] : elements_in_blocks_2[i]);
    }
    for(uint32_t i = 0; i < batch_size; i++){
      counter_indexes[i] = bucket_ids[i] % _buckets_per_block;
    }

    ar_counter bucket_start_indexes =
      exclusive_reduce_many<batch_size>(block_ids, counter_indexes);

    ar_counter counter_values;
    for(uint32_t i = 0; i < batch_size; i++){
      counter_values[i] = read_counter(block_ids[i],
        counter_indexes[i]);
      statuses[offset + i] = !((elements_in_blocks[i] == _max_fingerprints_per_block) ||
        (counter_values[i] == _slots_per_bucket));

      if(statuses[offset + i]){
        write_fingerprint_left_displace(_storage[block_ids[i]],
          bucket_start_indexes[i], fingerprints[i]);
        increment_fullness_counter(_storage[block_ids[i]],
          counter_indexes[i]);
      }
      // Adjust for stale count (-1)
      if(_block_fullness_array_enabled && (elements_in_blocks[i] ==
        _max_fingerprints_per_block - 1)){
        _block_fullness_array[block_ids[i]] = 1;
      }
    }
    // Set OTA if you placed the item in the secondary bucket/block
    for(uint32_t i = 0; i < batch_size; i++){
      if(_morton_filter_functionality_enabled && (statuses[offset + i]) &&
        (!try_first_block_insert[i])){
        set_overflow_status(bucket_ids_1[i], fingerprints[i],
          block_ids_1[i], counter_indexes[i]);
      }
    }
    // Resolve lingering collisions
    for(uint32_t i = 0; i < batch_size; i++){
      if(!statuses[offset + i]){
        statuses[offset + i] = random_kickout_cuckoo(bucket_ids_1[i],
          fingerprints[i]);
      }
    }
  }

  inline counter_t report_fsa_load(hash_t block_id){
    // WARNING: If the branch is taken, we're not actually finding the load,
    // rather whether the last slot of the fingerprint storage array
    // is empty or non-empty.  For first-fit this doesn't
    // matter because we only care if there is at least one free slot.
    // However, it does matter for two-choice.
    if(_special_null_fingerprint){
      bool is_nonzero_fingerprint = read_fingerprint(block_id,
        _max_fingerprints_per_block - 1);
      return is_nonzero_fingerprint ?
        _max_fingerprints_per_block : _max_fingerprints_per_block - 1;
    }
    else{
      return exclusive_reduce(_storage[block_id], _buckets_per_block);
    }
  }

  inline void first_level_store_many(const ar_hash& bucket_ids,
    const ar_atom& fingerprints, const ar_hash& block_ids,
    ar_counter& counter_indexes,
    ar_counter& bucket_start_indexes, ar_counter& elements_in_blocks,
    ar_counter& counter_values, std::vector<bool>& statuses, const hash_t offset){

    for(uint32_t i = 0; i < batch_size; i++){
      counter_indexes[i] = bucket_ids[i] % _buckets_per_block;
    }

    for(uint32_t i = 0; i < batch_size; i++){
      elements_in_blocks[i] = report_fsa_load(block_ids[i]);
    }
    bucket_start_indexes =
      exclusive_reduce_many<batch_size>(block_ids, counter_indexes);

    for(uint32_t i = 0; i < batch_size; i++){
      counter_values[i] = read_counter(block_ids[i], counter_indexes[i]);
      statuses[offset + i] = !((elements_in_blocks[i] ==
        _max_fingerprints_per_block) | (counter_values[i] == _slots_per_bucket));

      if(statuses[offset + i]){
        write_fingerprint_left_displace(_storage[block_ids[i]],
          bucket_start_indexes[i], fingerprints[i]);
        increment_fullness_counter(_storage[block_ids[i]],
          counter_indexes[i]);
      }
      // -1 because we just added something so the count is stale.
      if(_block_fullness_array_enabled && (elements_in_blocks[i] ==
        _max_fingerprints_per_block - 1)){
        _block_fullness_array[block_ids[i]] = 1;
      }
    }
  }

  template<uint64_t num_buckets>
  INLINE bool conflict_exists(BlockedBF::BloomFilter<num_buckets>& bf,
    const hash_t block_id){
    return bf.contains_and_update(block_id);
  }

  inline void table_store_many_two_choice(const ar_hash& bucket_ids_1,
    const ar_atom& fingerprints, std::vector<bool>& statuses,
    const hash_t offset){

    // Used if there is a conflict in the batch (two or more fingerprints that
    // would modify the same block)
    ar_store_params c1;  // Bucket/Block candidate 1
    ar_store_params c2;  // Bucket/Block candidate 2

    // These get used by two_choice_store_many
    ar_hash bucket_ids_2; // secondary
    ar_hash block_ids_1; // primary
    ar_hash block_ids_2; // secondary

    for(uint32_t i = 0; i < batch_size; i++){
      bucket_ids_2[i] = determine_alternate_bucket(bucket_ids_1[i],
        fingerprints[i]);
    }

    // Number of buckets in the Bloom filter
    constexpr uint64_t num_buckets = 128; // Must be power of 2
    BlockedBF::BloomFilter<num_buckets> bf;
    std::bitset<2 * batch_size> conflict_vector;
    for(uint32_t i = 0; i < batch_size; i++){
      block_ids_1[i] = bucket_ids_1[i] / _buckets_per_block;
      block_ids_2[i] = bucket_ids_2[i] / _buckets_per_block;
      conflict_vector[i << 1] = conflict_exists<num_buckets>(bf, block_ids_1[i]);
      conflict_vector[(i << 1) + 1] = conflict_exists<num_buckets>(bf, block_ids_2[i]);
    }

    if(_handle_conflicts && conflict_vector.any()){
      for(uint32_t i = 0; i < batch_size; i++){
        statuses[offset + i] = first_level_store(bucket_ids_1[i],
          fingerprints[i], c1[i]);
        if(!statuses[offset + i]){  // Try again
          statuses[offset + i] = first_level_store(bucket_ids_2[i],
            fingerprints[i], c2[i]);
          if(_morton_filter_functionality_enabled && statuses[offset + i]){
            set_overflow_status(bucket_ids_1[i], fingerprints[i],
              c1[i].block_id, c1[i].counter_index);
          }
        }
      }
    }
    else{ // Hopefully the common case
      return two_choice_store_many(bucket_ids_1, bucket_ids_2, fingerprints,
        block_ids_1, block_ids_2, statuses, offset);
    }

    for(uint_fast32_t i = 0; i < batch_size; i++){
      bool status1 = statuses[offset + i];
      InsertStatus status2 = InsertStatus::FAILED_TO_INSERT;

      if(_collision_resolution_enabled && !status1){
          status2 = resolve_collision(bucket_ids_1[i], bucket_ids_2[i],
            fingerprints[i], c1[i], c2[i]);
      }
      //Mark overflow here if item is placed in a secondary bucket
      if(_morton_filter_functionality_enabled &
        (status2 == InsertStatus::PLACED_IN_SECONDARY_BUCKET)){
        set_overflow_status(bucket_ids_1[i], fingerprints[i], block_ids_1[i],
          c1[i].counter_index);
      }

      statuses[offset + i] = status1 |
        (status2 != InsertStatus::FAILED_TO_INSERT);

      if(/*_DEBUG &&*/ !statuses[offset + i]){
      //std::cerr << "Table store failed on bucket " << bucket_ids[i] <<
       // " and fingerprint " << fingerprints[i] << ".\n" << std::endl;
      }
    }
  }

  inline void table_store_many(const ar_hash& bucket_ids,
    const ar_atom& fingerprints, std::vector<bool>& statuses,
    const hash_t offset){

    ar_store_params c1;  // Bucket/Block candidate 1
    ar_store_params c2;  // Bucket/Block candidate 2

    // These get passed in to first_level_store_many and are populated
    ar_hash block_ids;
    ar_counter counter_indexes;
    ar_counter bucket_start_indexes;
    ar_counter elements_in_blocks;
    ar_counter counter_values;

    // Number of buckets in the Bloom filter
    constexpr uint64_t num_buckets = 64; // Must be power of 2
    BlockedBF::BloomFilter<num_buckets> bf;
    std::bitset<batch_size> conflict_vector;
    for(uint32_t i = 0; i < batch_size; i++){
      block_ids[i] = bucket_ids[i] / _buckets_per_block;
      conflict_vector[i] = conflict_exists<num_buckets>(bf, block_ids[i]);
    }

    if(__builtin_expect(_handle_conflicts && conflict_vector.any(), 0)){
      for(uint32_t i = 0; i < batch_size; i++){
        statuses[offset + i] = first_level_store(bucket_ids[i],
          fingerprints[i], c1[i]);
      }
    }
    else{
      first_level_store_many(bucket_ids, fingerprints, block_ids,
        counter_indexes, bucket_start_indexes, elements_in_blocks,
        counter_values, statuses, offset);
      for(uint32_t i = 0; i < batch_size; i++){
        c1[i].block_id = block_ids[i];
        c1[i].counter_index = counter_indexes[i];
      }
    }
    for(uint_fast32_t i = 0; i < batch_size; i++){
      bool status1 = statuses[offset + i];
      bool status2 = false;
      InsertStatus status3 = InsertStatus::FAILED_TO_INSERT;

      if(_remap_enabled && !status1){
        hash_t secondary_bucket_id = determine_alternate_bucket(bucket_ids[i],
          fingerprints[i]);
        status2 = first_level_store(secondary_bucket_id, fingerprints[i], c2[i]);

        if(_collision_resolution_enabled && !status2){
          c1[i].counter_value = read_counter(c1[i].block_id,
          c1[i].counter_index);
          c1[i].elements_in_block = report_fsa_load(c1[i].block_id);
          c1[i].bucket_start_index = get_bucket_start_index(c1[i].block_id,
            c1[i].counter_index);
          status3 = resolve_collision(bucket_ids[i], secondary_bucket_id,
            fingerprints[i], c1[i], c2[i]);
        }
      }
      //Mark overflow here if item is placed in a secondary bucket
      if(_morton_filter_functionality_enabled &
        (status2 | (status3 == InsertStatus::PLACED_IN_SECONDARY_BUCKET))){
        set_overflow_status(bucket_ids[i], fingerprints[i], c1[i].block_id,
          c1[i].counter_index);
      }

      statuses[offset + i] = status1 | status2 |
        (status3 != InsertStatus::FAILED_TO_INSERT);

      if(/*_DEBUG &&*/ !statuses[offset + i]){
      //std::cerr << "Table store failed on bucket " << bucket_ids[i] <<
       // " and fingerprint " << fingerprints[i] << ".\n" << std::endl;
      }
    }
  }

  inline bool first_level_store(hash_t bucket_id, atom_t fingerprint,
    StoreParams& sp){
    /*if(_block_fullness_array_enabled && _block_fullness_array[bucket_id]){
      return false;
    } This won't work here */

    sp.block_id = bucket_id / _buckets_per_block;
    sp.counter_index = (bucket_id % _buckets_per_block);

    switch(_reduction_method){
      case ReductionMethodEnum::NAIVE_FULL_EXCLUSIVE_SCAN:{
        counter_t* exclusive_scan_ptr = full_exclusive_scan(sp.block_id);
        sp.bucket_start_index = exclusive_scan_ptr[sp.counter_index];
        sp.elements_in_block = exclusive_scan_ptr[_buckets_per_block];
      }
        break;

      case ReductionMethodEnum::PARALLEL_REDUCE:
        sp.bucket_start_index = exclusive_reduce_with_parallel_sum(
          _storage[sp.block_id], sp.counter_index);
        sp.elements_in_block = report_fsa_load(sp.block_id);
        break;

      case ReductionMethodEnum::POP_CNT:
        sp.bucket_start_index = exclusive_reduce_with_popcount(
          _storage[sp.block_id], sp.counter_index);
        sp.elements_in_block = report_fsa_load(sp.block_id);
        break;
      default:
        std::cerr << "ERROR: Undefined setting for _reduction_method\n";
        exit(1);
        break;
    }

    if(_DEBUG){
      std::cout << "Inserting fingerprint " << fingerprint << " into block "
        << sp.block_id << ", which has " <<
        static_cast<uint32_t>(sp.elements_in_block) <<
        " elements in the block\n";
    }

    sp.counter_value = read_counter(sp.block_id, sp.counter_index);

    bool remap_necessary = (sp.elements_in_block == _max_fingerprints_per_block)
      | (sp.counter_value == _slots_per_bucket);

    if(remap_necessary) return false;

    if(_DEBUG){
      std::cout << "Inserting fingerprint " << fingerprint << " at offset " <<
        _fingerprint_offset << " bits from start of bucket at index " <<
        static_cast<uint32_t>(sp.bucket_start_index) << ". Fingerprint is "
        << _fingerprint_len_bits <<
        " wide. " << std::endl;
    }

    // Shift everything over by one entry starting at the start of the bucket
    // Stuff to the right of adjusted_index remains unmodified
    write_fingerprint_left_displace(_storage[sp.block_id], sp.bucket_start_index,
      fingerprint);

    // Increment the counter
    increment_fullness_counter(_storage[sp.block_id], sp.counter_index);

    // Update block fullness array (-1 because value is stale)
    if(_block_fullness_array_enabled && (sp.elements_in_block ==
      _max_fingerprints_per_block - 1)){
      _block_fullness_array[sp.block_id] = true;
    }
    return true;
  }


  // Relocation on bucket overflows
  inline bool try_relocation(hash_t bucket_id, atom_t fingerprint,
    const StoreParams& sp){
    for(counter_t slot_id = 0; slot_id < sp.counter_value; slot_id++){
        atom_t candidate_fingerprint_to_evict =
          read_fingerprint(sp.block_id, sp.bucket_start_index + slot_id);
        hash_t alternate_bucket_id = determine_alternate_bucket(bucket_id,
          candidate_fingerprint_to_evict);

        bool proceed_to_first_level_store = true;
        if(_block_fullness_array_enabled){
          proceed_to_first_level_store = ! _block_fullness_array[alternate_bucket_id / _buckets_per_block];
        }

        StoreParams sp_prime;
        // Check if we can remap the candidate fingerprint in its alternate
        // candidate bucket
        if(proceed_to_first_level_store &&
          first_level_store(alternate_bucket_id,
          candidate_fingerprint_to_evict, sp_prime)){
          // Write fingerprint into newly vacated slot
          write_fingerprint(_storage[sp.block_id], sp.bucket_start_index +
            slot_id, fingerprint);

          // Mark the overflow
          if(_morton_filter_functionality_enabled){
            set_overflow_status(bucket_id, candidate_fingerprint_to_evict,
              sp.block_id, sp.counter_index);
          }

          return true;
        }
    }
    return false;
  }

  // This function kicks an item out of another bucket to make room for an
  // item at index bucket id.  Both items are within the same block.
  // This function with some minor modifications may also be conservatively
  // applied to bucket overflows. To do so, it needs to change the range of
  // values that get checked to a single bucket rather than a block.

  // However, this function does more work because it requires performing a
  // deletion plus an add with shifting in the initial block rather than just
  // a simple eviction and replacement.
  INLINE bool try_relocation_on_block_overflow(hash_t bucket_id,
    atom_t fingerprint, StoreParams& sp){ //* Need to pass in loop parameters
    // Get the first bucket in the block
    bucket_id -= (bucket_id % _buckets_per_block); //* Do before calling fcn
    for(counter_t counter_index = 0; counter_index < _buckets_per_block;
      counter_index++, bucket_id++){ //* Bucket overflow, this loop has 1 iter.
      // Get the start of the block
      counter_t bucket_start_index = exclusive_reduce(
        _storage[sp.block_id], counter_index);
      counter_t occupied_slots = read_counter(sp.block_id, counter_index);
      for(counter_t slot_id = 0; slot_id < occupied_slots; slot_id++){
        atom_t candidate_fingerprint_to_evict =
          read_fingerprint(sp.block_id, bucket_start_index + slot_id);
        hash_t alternate_bucket_id = determine_alternate_bucket(bucket_id,
          candidate_fingerprint_to_evict);
        StoreParams sp_prime;
        bool proceed_to_first_level_store = true;
        if(_block_fullness_array_enabled){
          proceed_to_first_level_store = ! _block_fullness_array[alternate_bucket_id / _buckets_per_block];
        }
        // Check if we can remap the candidate fingerprint to its alternate
        // candidate bucket
        if(proceed_to_first_level_store &&
          first_level_store(alternate_bucket_id,
          candidate_fingerprint_to_evict, sp_prime)){
          // Delete element that was remapped (i.e., don't store it twice).
          delete_fingerprint_right_displace(_storage[sp.block_id],
            bucket_start_index + slot_id);
          decrement_fullness_counter(_storage[sp.block_id], counter_index,
            occupied_slots);
          // Recompute the bucket start index for the fingerprint we want to
          // insert.  It will have be one less if we deleted an item in a
          // bucket that precedes the bucket where we want to make an insertion.
          // Otherwise, there is no change.
          sp.bucket_start_index = (counter_index < sp.counter_index ?
            sp.bucket_start_index - 1: sp.bucket_start_index);

          write_fingerprint_left_displace(_storage[sp.block_id],
            sp.bucket_start_index, fingerprint);
          increment_fullness_counter(_storage[sp.block_id], sp.counter_index);

          // Mark the overflow bit
          if(_morton_filter_functionality_enabled){
            set_overflow_status(bucket_id, candidate_fingerprint_to_evict,
              sp.block_id, sp.counter_index);
          }

          if(_DEBUG) std::cout << "Exiting block overflow with success.\n";
          return true;
        }

      }
    }
    if(_DEBUG) std::cout << "Exiting block overflow with failure.\n";
    return false;
  }

  INLINE hash_t get_ota_index(hash_t bucket_id, atom_t fingerprint) const{
    constexpr bool scale_bucket_id = _resizing_enabled &
      (_morton_ota_hashing_method !=
      OverflowTrackingArrayHashingMethodEnum::LEMIRE_FINGERPRINT_MULTIPLY) &
      (__builtin_popcountll(_buckets_per_block) != 1);
    if(scale_bucket_id){
      hash_t old_block_id = (bucket_id / _buckets_per_block) >> _resize_count;
      hash_t lbi = old_block_id % _buckets_per_block;
      bucket_id = old_block_id * _buckets_per_block + lbi;
    }
    hash_t ota_index;
    switch(_morton_ota_hashing_method){
      case OverflowTrackingArrayHashingMethodEnum::LEMIRE_FINGERPRINT_MULTIPLY:
        ota_index = util::fast_mod_alternative<atom_t>(fingerprint,
         _ota_len_bits, _fingerprint_len_bits);
        break;

      case OverflowTrackingArrayHashingMethodEnum::RAW_BUCKET_HASH:
        ota_index = bucket_id % _ota_len_bits;
        break;

      case OverflowTrackingArrayHashingMethodEnum::CLUSTERED_BUCKET_HASH:{
        // Correct for compile-time division by zero when this isn't called
        // because the OTA's length is zero bits.
        constexpr uint64_t divisor = _ota_len_bits > 0 ?
          (_buckets_per_block + _ota_len_bits - 1) /_ota_len_bits :
          1;
        ota_index = (bucket_id % _buckets_per_block) / divisor;
        break;
      }
    }
    return ota_index;
  }


  // This is designed to work with OTAs that are powers of two bits in length.
  inline bool check_bloom_filter_ota(hash_t bucket_id, atom_t fingerprint,
    const hash_t block_id) const{
    const uint_fast64_t num_hash_fcns = 2;
    std::bitset<num_hash_fcns> status;
    const atom_t hash = bucket_id * 0xff51afd7ed558ccdULL; // MurmurHash mixing
    for(uint64_t i = 0; i < num_hash_fcns; i++){
      // Pass in hash rather than true bucket ID
      hash_t ota_index = get_ota_index(hash >> (i * _ota_len_bits),
        fingerprint);
      status[i] = _storage[block_id].read_bit(_overflow_tracking_array_offset +
        ota_index);
    }
    return status.all();
  }

  // This is designed to work with OTAs that are powers of two bits in length.
  inline void set_bloom_filter_ota(const hash_t bucket_id, atom_t fingerprint,
    const hash_t block_id){
    constexpr uint_fast64_t num_hash_fcns = 2;
    const atom_t hash = bucket_id * 0xff51afd7ed558ccdULL; // MurmurHash mixing
    for(uint64_t i = 0; i < num_hash_fcns; i++){
      // Pass in hash rather than true bucket ID
      hash_t ota_index = get_ota_index(hash >> (i * _ota_len_bits), fingerprint);
      _storage[block_id].sticky_set_bit(_overflow_tracking_array_offset +
        ota_index, 1);
    }
  }

  // Check whether an OTA bit is set
  inline bool get_overflow_status(hash_t bucket_id, atom_t fingerprint) const{
    hash_t ota_index = get_ota_index(bucket_id, fingerprint);
    hash_t block_id = bucket_id / _buckets_per_block;
    // Bloom filter
    if(_use_bloom_ota){ // Not yet implemented for selective Morton filter
      return check_bloom_filter_ota(bucket_id, fingerprint, block_id);
    }
    // Bit vector
    constexpr bool using_selective_mf = _ota_lbi_insertion_threshold > - 1;
    if(using_selective_mf){
      hash_t lbi = bucket_id % _buckets_per_block;
      return _storage[block_id].read_bit(_overflow_tracking_array_offset +
        ota_index) | (static_cast<int16_t>(lbi) <= _ota_lbi_insertion_threshold);
    }
    else{
      return _storage[block_id].read_bit(_overflow_tracking_array_offset +
        ota_index);
    }
  }

  // Set OTA bit on overflow if not already set
  inline void set_overflow_status(const hash_t bucket_id,
    const atom_t fingerprint, const hash_t block_id, const hash_t){
    // Bloom filter
    if(_use_bloom_ota){ // Not yet implemented for selective Morton filter
      return set_bloom_filter_ota(bucket_id, fingerprint, block_id);
    }
    // Bit vector
    hash_t ota_index = get_ota_index(bucket_id, fingerprint);
    constexpr bool using_selective_mf = _ota_lbi_insertion_threshold > -1;
    if(using_selective_mf){
      _storage[block_id].sticky_set_bit(
        _overflow_tracking_array_offset + ota_index,
        static_cast<int16_t>(bucket_id % _buckets_per_block) > _ota_lbi_insertion_threshold);
    }
    else{
      _storage[block_id].sticky_set_bit(
        _overflow_tracking_array_offset + ota_index, 1);
    }
  }

  inline bool random_kickout_cuckoo(hash_t bucket_id, atom_t fingerprint){
    uint_fast16_t max_count = 300; // 40 works for many configurations
    uint_fast16_t count = 1;

    AccessCounter access_count;
    while(count <= max_count){
      // Compute block ID from candidate bucket
      uint64_t block_id = bucket_id / _buckets_per_block;
      // Compute block-local bucket ID lbi
      const uint16_t counter_index = (bucket_id % _buckets_per_block);
      // Use the block-local bucket ID to read the fullness counter at FCA[lbi]
      counter_t full_slots = read_counter(block_id, counter_index);
      // Compute offset into fingerprint storage array using a summed reduction
      counter_t bucket_start_index = get_bucket_start_index(block_id,
        counter_index);
      // Compute whether the block's FSA is full with another reduction
      counter_t fsa_load_in_fingerprints = get_bucket_start_index(block_id,
        _buckets_per_block);

      // If FCA[lbi] < MAX_SLOTS_PER_BUCKET and the total load < MAX_FINGERPRINTS
      if(__builtin_expect((full_slots < _slots_per_bucket) &
        (fsa_load_in_fingerprints < _max_fingerprints_per_block), 0)){
        // Perform a first level store of the fingerprint into the bucket
        // 1) Add the fingerprint to the bucket within the block
        write_fingerprint_left_displace(_storage[block_id], bucket_start_index,
          fingerprint);
        // 2) Increment the fullness counter associated with the bucket
        increment_fullness_counter(_storage[block_id], counter_index);
        // 3) Update the block fullness array (-1 because value is stale)
        if(_block_fullness_array_enabled &&
          fsa_load_in_fingerprints == _max_fingerprints_per_block - 1){
          _block_fullness_array[block_id] = true;
        }
        if(_print_access_counts){
          access_count.primary_count++;
          std::cout << (uint32_t) access_count.primary_count + 1 << ","
            << (uint32_t) access_count.block_overflow_count << ","
            << (uint32_t) access_count.bucket_overflow_count << ","
            << (uint32_t) access_count.block_and_bucket_overflow_count <<
            std::endl;
        }
        return true;
      }
      // Bucket Overflow (takes precedence over block overflows)
      // Else If the fullness counter at FCA[lbi] is has a value of max slots
      // per bucket
      else if(__builtin_expect(full_slots == _slots_per_bucket, 0)){
        // undergo collision resolution for bucket overflows
        // 1) Select a random fingerprint from the bucket to evict and save it
        // as f2 off to the side
        uint32_t slot_id = _slots_per_bucket - 1; //rand() % _slots_per_bucket;
        atom_t f2 = read_fingerprint(block_id, bucket_start_index + slot_id);
        // 1b) Set the associated OTA bit that corresponds to jettisoning f2
        if(_morton_filter_functionality_enabled){
          set_overflow_status(bucket_id, f2, block_id, counter_index);
        }
        // 2) Overwrite f2 in the bucket with the new fingerprint "fingerprint"
        write_fingerprint(_storage[block_id], bucket_start_index + slot_id,
          fingerprint);
        // 3) Compute the alternate bucket for f2 and set bucket_id to this
        // bucket
        bucket_id = determine_alternate_bucket(bucket_id, f2);
        // 4) Set fingerprint to f2
        fingerprint = f2;
        if(_print_access_counts){
          if(fsa_load_in_fingerprints == _max_fingerprints_per_block){
            access_count.block_and_bucket_overflow_count++;
          }
          else{
            access_count.bucket_overflow_count++;
          }
        }
      }

      // Block Overflow (Only come here if you are not also a bucket overflow)
      // Else the block's FSA is full, then undergo collision resolution for
      // block overflows
      else{
        // 1) Select a random fingerprint f3 from the block to cache off to the
        // side
        // 1a) First select a bucket that has elements in it
        hash_t eviction_bucket_id; // Block local bucket ID not global
        counter_t eviction_bucket_count;
        // Slow but functional way to get the block-local bucket ID of a
        // non-empty bucket
        do{
          eviction_bucket_id = rand() % (_buckets_per_block);
          eviction_bucket_count = read_counter(block_id, eviction_bucket_id);
        } while(eviction_bucket_count == 0);
        // 1b) Select a random slot and cache the fingerprint f3
        hash_t eviction_slot_id = rand() % eviction_bucket_count;
        hash_t eviction_bucket_offset = get_bucket_start_index(block_id,
          eviction_bucket_id);
        atom_t f3 = read_fingerprint(block_id,
          eviction_bucket_offset + eviction_slot_id);
        // 2) Delete f3 from the block
        delete_fingerprint_right_displace(_storage[block_id],
          eviction_bucket_offset + eviction_slot_id);
        // 2b) Set the associated OTA bit for f3
        if(_morton_filter_functionality_enabled){
          hash_t full_eviction_bucket_id = block_id * _buckets_per_block +
            eviction_bucket_id;
          set_overflow_status(full_eviction_bucket_id, f3, block_id,
            eviction_bucket_id);
        }
        // 3) Decrement the fullness counter associated with f3 if delete
        // doesn't already
        decrement_fullness_counter(_storage[block_id], eviction_bucket_id,
          eviction_bucket_count);
        // 4) Insert fingerprint 'fingerprint' into its bucket and increment
        // its fullness counter
        // 4a) Correct initial bucket offset if deleting an item would've
        // messed it up
        bucket_start_index = (eviction_bucket_id < counter_index ?
          bucket_start_index - 1 : bucket_start_index);
        // 4b) Insert
        write_fingerprint_left_displace(_storage[block_id], bucket_start_index,
          fingerprint);
        // 4c) Fullness counter increment
        increment_fullness_counter(_storage[block_id], counter_index);
        // 5) Compute the alternate bucket for f3 and set bucket_id to
        // this bucket
        bucket_id = determine_alternate_bucket(block_id * _buckets_per_block +
          eviction_bucket_id, f3);
        // 6) Set fingerprint to f3
        fingerprint = f3;
        if(_print_access_counts) access_count.block_overflow_count++;
      } // End of block overflow resolution
      count ++;
    } // End of while loop
    // If you exit the while loop here, it means that the max count has been
    // exceeded.
    return false;
  }

  inline void double_capacity(){
    resize<1>();
  }

  inline void quadruple_capacity(){
    resize<2>();
  }

  inline void octuple_capacity(){
    resize<3>();
  }

  // Resizes the filter by increasing the number of blocks in the filter by factor of pow(2, log2_resize).
  // This code assumes the number of buckets is even and that there are no partial blocks.  This assumption
  // is made elsewhere in the code and is enforced by the constructor.
  template<uint64_t log2_resize>
  inline void resize(){
    if(!_resizing_enabled){
      std::cerr << "Set the _resizing_enabled flag to use the resize() or double_capacity() methods\n";
      exit(1);
    }

    constexpr hash_t one = 1;
    constexpr uint64_t resize_factor = (one << log2_resize);
    uint64_t new_total_buckets = _total_buckets * resize_factor;
    uint64_t new_total_slots = new_total_buckets * _slots_per_bucket; // Virtual not actual
    uint64_t new_total_blocks = (new_total_slots + _slots_per_bucket * _buckets_per_block - 1) / (_slots_per_bucket * _buckets_per_block); // Round up to next whole block
    //std::vector<bool> new_block_fullness_array(_block_fullness_array_enabled ? new_total_blocks : 0, 0);
    block_t* new_storage = allocate_cache_aligned_storage(new_total_blocks);
    block_t* old_storage = _storage;

    for(uint64_t block_id = 0; block_id < _total_blocks; block_id++){
      // Tracks how many fingerprints have been added already to each of the
      // two child blocks (0 at the start) but up to _fingerprints_per_block
      hash_t fsa_pointers[resize_factor]{};
      hash_t fsa_index = 0;
      // counter_index is synonymous with lbi from the VLDB'18 paper
      for(uint64_t counter_index = 0; counter_index < _buckets_per_block;
        counter_index++){
        counter_t fullness_counter = read_counter(block_id, counter_index);
        for(uint64_t slot_id = 0; slot_id < fullness_counter; slot_id++, fsa_index++){
          atom_t fingerprint = read_fingerprint(block_id, fsa_index);
          hash_t fsa_pointer_index = ((fingerprint >> (_fingerprint_len_bits - _resize_count - log2_resize)) & (resize_factor - one));
          hash_t new_block = resize_factor * block_id + fsa_pointer_index;
          write_fingerprint(new_storage[new_block],
            fsa_pointers[fsa_pointer_index], fingerprint);
          fsa_pointers[fsa_pointer_index]++;
          increment_fullness_counter(new_storage[new_block], counter_index);
        }
      }
      // Copy OTA to new blocks
      atom_t ota = old_storage[block_id].read_cross(_overflow_tracking_array_offset, _ota_len_bits, 0);
      uint64_t base_id = block_id * resize_factor;
      for(uint64_t new_block_id = base_id; new_block_id < base_id + resize_factor; new_block_id++){
        new_storage[new_block_id].add_cross(_overflow_tracking_array_offset, _ota_len_bits, 0, ota);
      }
    }
    _total_buckets = new_total_buckets;
    _total_slots = new_total_slots;
    _total_blocks = new_total_blocks;
    _storage = new_storage;
    // TODO: Finish implementing _block_fullness_array = new_block_fullness_array;
    _resize_count+=log2_resize;
    free(old_storage); // FIXME: Only works with the aligned_alloc call
  }

  // The main function for resolving collisions during insertions.  It does
  // a two level breadth-first search but then reverts to using a Morton-
  // filter-specific variant of Fan et al.'s random kickout
  // algorithm (random_kickout_cuckoo).
  // I did this to avoid storing explored paths and keeping track
  // of visited nodes.  It didn't seem worth the overhead since the
  // vast majority of collisions can be resolved by looking solely at
  // the primary and secondary blocks.  For those that remain, most
  // are resolvable by going just one level
  // deeper (try_relocation_on_block_overflow and
  // try_relocation_on_bucket_overflow).
  INLINE InsertStatus resolve_collision(hash_t bucket_id,
    hash_t secondary_bucket_id,
    atom_t fingerprint, StoreParams& c1, StoreParams& c2){

    // Lookups are faster when you leave this commented out, since
    // it improves biasing in favor of the primary hash function.
    // However, random_kickout_cuckoo is a stand-alone insertion algorithm.
    //return random_kickout_cuckoo(bucket_id, fingerprint) ? InsertStatus::PLACED_IN_PRIMARY_BUCKET : InsertStatus::FAILED_TO_INSERT;

    // Pure block overflow?
    bool primary_block_overflow = c1.counter_value != _slots_per_bucket;
    bool secondary_block_overflow = c2.counter_value != _slots_per_bucket;

    // Try this first since most overflows are block overflows
    if(primary_block_overflow &&
      try_relocation_on_block_overflow(bucket_id, fingerprint, c1)){
      if(_DEBUG) std::cout << "Success on third relocation attempt\n";
      return InsertStatus::PLACED_IN_PRIMARY_BUCKET;
    }

    // Bucket overflow
    if(/*(!primary_block_overflow) &&*/ try_relocation(bucket_id, fingerprint, c1)){
      if(_DEBUG) std::cout << "Success on first relocation attempt\n";
      return InsertStatus::PLACED_IN_PRIMARY_BUCKET;
    }

    // Pure block overflow?
    if(secondary_block_overflow &&
      try_relocation_on_block_overflow(secondary_bucket_id, fingerprint, c2)){
      if(_DEBUG) std::cout << "Success on fourth relocation attempt\n";
      return InsertStatus::PLACED_IN_SECONDARY_BUCKET;
    }

    if(/*(!secondary_block_overflow) &&*/
      try_relocation(secondary_bucket_id, fingerprint, c2)){
      if(_DEBUG) std::cout << "Success on second relocation attempt\n";
      return InsertStatus::PLACED_IN_SECONDARY_BUCKET;
    }

    InsertStatus return_status = random_kickout_cuckoo(bucket_id, fingerprint) ?
      InsertStatus::PLACED_IN_PRIMARY_BUCKET :
      InsertStatus::FAILED_TO_INSERT;
    return return_status;

    return InsertStatus::FAILED_TO_INSERT;
  }

  // Store the fingerprint in the bucket specified by bucket_id
  inline bool table_store(hash_t bucket_id, atom_t fingerprint){
    StoreParams c1;  // Bucket/Block candidate 1
    StoreParams c2;  // Bucket/Block candidate 2
    bool status1 = first_level_store(bucket_id, fingerprint, c1);
    bool status2 = false;
    InsertStatus status3 = InsertStatus::FAILED_TO_INSERT;
    if(_remap_enabled && !status1){
      hash_t secondary_bucket_id = determine_alternate_bucket(bucket_id,
        fingerprint);
      status2 = first_level_store(secondary_bucket_id, fingerprint, c2);
      if(_collision_resolution_enabled && !status2){
        status3 = resolve_collision(bucket_id, secondary_bucket_id,
          fingerprint, c1, c2);
      }
    }
    //Mark overflow here if item is placed in a secondary bucket
    if(_morton_filter_functionality_enabled &
      (status2 | (status3 == InsertStatus::PLACED_IN_SECONDARY_BUCKET))){
      set_overflow_status(bucket_id, fingerprint, c1.block_id, c1.counter_index);
    }

    bool net_status = status1 | status2 |
      (status3 != InsertStatus::FAILED_TO_INSERT);

    if(_print_access_counts && (status1 | status2)){
      std::cout << status1 + (status2 << 1) << ",0,0,0" << std::endl;
    }

    if(/*_DEBUG &&*/ !net_status){
      std::cerr << "Table store failed on bucket " << bucket_id <<
        " and fingerprint " << fingerprint << ".\n" << std::endl;
    }
    return net_status;
  }

  // Returns the layout of the table as it is actually stored
  std::string physical_layout(){
    std::stringstream ss;
    for(hash_t block_id = 0; block_id < _total_blocks; block_id++){
      ss << "Block " << block_id << ": [ ";
      for(int32_t slot = _max_fingerprints_per_block - 1; slot > -1; slot--){
        ss << read_fingerprint(block_id, slot);
      }
      ss << "]" << std::endl << "[ ";
      for(int32_t counter = _buckets_per_block - 1; counter > - 1; counter--){
        ss << read_counter(block_id, counter);
      }
      ss << "]" << std::endl;
    }
    return ss.str();
  }

  std::string as_string(){
    std::stringstream ss;
    hash_t buckets_processed = 0;
    for(hash_t b_id = 0; b_id < _total_blocks; b_id++){
      // Prefix sum the block's counters
      counter_t* scan = full_exclusive_scan(b_id);
      // Complex check is necessary if we round up to the next full block
      // but don't use all of the buckets. This might be mandated by the
      // hashing scheme that we use.
      for(hash_t bucket = 0; bucket < _buckets_per_block
        && buckets_processed < _total_buckets; bucket++){
      //for(uint64_t bucket = 0; bucket < cf._buckets_per_block; bucket++){
        ss << "[ ";
        // Figure out where to start reading in the block
        uint16_t read_index = scan[bucket];
        uint16_t occupied_slots = read_counter(b_id, bucket);
        for(uint16_t slot = 0; slot < _slots_per_bucket; slot++){
          // Reading an item at a time is inefficient, but this code doesn't
          // need to be fast.  The control flow isn't elegant either.
          atom_t fingerprint = read_fingerprint(b_id, read_index + slot);
          if(slot < occupied_slots){
            ss << fingerprint << " ";
          }
          else{
            ss << "EMPTY ";
          }
        }
        ss << "]" << " " << occupied_slots << std::endl;
        buckets_processed++;
       }
     }
     return ss.str();
  }

  // Query the actual max slot compression ratio
  // This figure is the quotient of the capacity in fingerprints per block
  // and the net logical slots within a block.
  double get_max_slot_compression_ratio(){
    return static_cast<double>(_max_fingerprints_per_block)/(_slots_per_bucket *
      _buckets_per_block);
  }

  void print_bucket_and_block_load_histograms(){
    uint64_t bucket_load_histogram[_slots_per_bucket + 1]{}; // 0 to S (see VLDB)
    uint64_t block_load_histogram[_max_fingerprints_per_block + 1]{}; // 0 to B
    uint64_t combined_histogram[_slots_per_bucket + 1][_max_fingerprints_per_block + 1]{};
    for(uint64_t b_id = 0; b_id < _total_blocks; b_id++){
      uint64_t block_load = get_bucket_start_index(b_id, _buckets_per_block);
        block_load_histogram[block_load]++;
      for(uint64_t bucket_id = 0; bucket_id < _buckets_per_block; bucket_id++){
        uint64_t bucket_load = read_counter(b_id, bucket_id);
        bucket_load_histogram[bucket_load]++;
        combined_histogram[bucket_load][block_load]++;
      }
    }
    uint64_t buckets_with_spare_capacity = 0;
    // Print bucket load histogram
    std::cout << "TABLE";
    for(uint32_t bucket_load = 0; bucket_load <= _slots_per_bucket;
      bucket_load++){
      std::cout << "\t" << bucket_load;
    }
    std::cout << std::endl;
    std::cout << "BUCKET_LOAD_HISTOGRAM";
    for(uint32_t bucket_load = 0; bucket_load <= _slots_per_bucket;
      bucket_load++){
      std::cout << "\t" << bucket_load_histogram[bucket_load];
      if(bucket_load < _slots_per_bucket){
        buckets_with_spare_capacity += bucket_load_histogram[bucket_load];
      }
    }
    std::cout << std::endl << std::endl;

    uint64_t blocks_with_spare_capacity = 0;
    // Print block load histogram
    std::cout << "TABLE";
    for(uint32_t block_load = 0; block_load <= _max_fingerprints_per_block;
      block_load++){
      std::cout << "\t" << block_load;
      if(block_load < _max_fingerprints_per_block){
        blocks_with_spare_capacity += block_load_histogram[block_load];
      }
    }
    std::cout << std::endl;
    std::cout << "BLOCK_LOAD_HISTOGRAM";
    for(uint32_t block_load = 0; block_load <= _max_fingerprints_per_block;
      block_load++){
      std::cout << "\t" << block_load_histogram[block_load];
    }
    std::cout << std::endl << std::endl;

    // Print combined load histogram
    std::cout << "BLOCK_AND_BUCKET_LOAD_HISTOGRAM\n";
    std::cout << "Fingerprints";
    for(uint32_t bucket_load = 0; bucket_load <= _slots_per_bucket;
      bucket_load++){
      std::cout << "\t" << bucket_load;
    }
    std::cout << std::endl;
    uint64_t buckets_with_spare_capacity_in_nonfull_block = 0;
    for(uint32_t block_load = 0; block_load <= _max_fingerprints_per_block;
      block_load++){
      std::cout << block_load;
      for(uint32_t bucket_load = 0; bucket_load <= _slots_per_bucket;
        bucket_load++){
        std::cout << "\t" << combined_histogram[bucket_load][block_load];
        if(block_load != _max_fingerprints_per_block &&
          bucket_load != _slots_per_bucket){
          buckets_with_spare_capacity_in_nonfull_block +=
            combined_histogram[bucket_load][block_load];
        }
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Percent of blocks with spare capacity: " <<
      100 * blocks_with_spare_capacity / static_cast<double>(_total_blocks) <<
      std::endl;

    std::cout << "Percent of buckets with spare capacity: " <<
      100 * buckets_with_spare_capacity / static_cast<double>(_total_buckets) <<
      std::endl;

    std::cout << "Percent of buckets with spare capacity in blocks with spare "
      "capacity: " << 100 * buckets_with_spare_capacity_in_nonfull_block /
      static_cast<double>(_total_buckets) << std::endl;
  }

  // I used this for some debugging that I was doing.  It supposed to match the output from
  // another program written in Python.
  std::string get_array_dimensions_as_string(){
    std::stringstream ss;
    ss << "OTA: " << _ota_len_bits << " x 1 bits = " << _ota_len_bits << "\n";
    ss << "FCA: " << _buckets_per_block << " x " << _fullness_counter_width <<
      " bits = " << _buckets_per_block * _fullness_counter_width << "\n";
    ss << "FSA: " << _max_fingerprints_per_block << " x " <<
      _fingerprint_len_bits << " bits = " <<
      _max_fingerprints_per_block * _fingerprint_len_bits << "\n";
    ss << "Net bits: " << _ota_len_bits << " + " << _buckets_per_block *
      _fullness_counter_width << " + " << _max_fingerprints_per_block *
      _fingerprint_len_bits  << " = " << _ota_len_bits + _buckets_per_block *
      _fullness_counter_width + _max_fingerprints_per_block *
      _fingerprint_len_bits << "\n";
    ss << "C: " << get_max_slot_compression_ratio() << "\n";
    return ss.str();
  }

  friend std::ostream& operator<<(std::ostream& os,
    CompressedCuckooFilter& cf){
    os << "CompressedCuckooFilter at address: " << &cf << std::endl;
    os << "Fingerprint length in bits: " << /*cf.*/cf._fingerprint_len_bits << std::endl;
    os << "Overflow tracking array length in bits: " << cf._ota_len_bits << std::endl;
    os << "Block length in bits: " << cf._block_size_bits << std::endl;
    os << "Virtual slots per bucket: " << cf._slots_per_bucket << std::endl;
    os << "Max slot compression ratio: " << cf.get_max_slot_compression_ratio()
      << std::endl;
    os << "Total slots: " << cf._total_slots << std::endl;
    os << "Total buckets: " << cf._total_buckets << std::endl;
    os << "Buckets per block: " << cf._buckets_per_block << std::endl;
    os << "Total blocks: " << cf._total_blocks << std::endl;
    os << "Max Fingerprints per block: " << cf._max_fingerprints_per_block << std::endl;
    os << "Fullness counters offset: " << cf._fullness_counters_offset << std::endl;
    os << "Reserved fields offset: " << cf._overflow_tracking_array_offset << std::endl;
    os << "Fingerprint offset: " << cf._fingerprint_offset << std::endl;

    if(_DEBUG){os << cf.as_string();}
    return os;
  }
};

} // End of CompressedCuckoo namespace

#endif
