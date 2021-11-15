#ifndef _VQF_FILTER_C_
#define _VQF_FILTER_C_

// https://github.com/splatlab/vqf

/*
 * ============================================================================
 *
 *       Filename:  vqf_filter.c
 *
 *         Author:  Prashant Pandey (), ppandey@berkeley.edu
 *   Organization:  LBNL/UCB
 *
 * ============================================================================
 */

#include <algorithm>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>  // portable to all x86 compilers
#include <tmmintrin.h>

#include "vqf_filter.h"
#include "vqf_precompute.h"

namespace vqfilter {

// ALT block check is set of 75% of the number of slots
#if TAG_BITS == 8
#define TAG_MASK 0xff
#define QUQU_SLOTS_PER_BLOCK 48
#define QUQU_BUCKETS_PER_BLOCK 80
#define QUQU_CHECK_ALT 92
#elif TAG_BITS == 12
#define TAG_MASK 0xfff
#define QUQU_SLOTS_PER_BLOCK 32
#define QUQU_BUCKETS_PER_BLOCK 96
#define QUQU_CHECK_ALT 104
#elif TAG_BITS == 16
#define TAG_MASK 0xffff
#define QUQU_SLOTS_PER_BLOCK 28 
#define QUQU_BUCKETS_PER_BLOCK 36
#define QUQU_CHECK_ALT 43 
#endif

#ifdef __AVX512BW__
extern __m512i SHUFFLE [];
extern __m512i SHUFFLE_REMOVE [];
extern __m512i SHUFFLE16 [];
extern __m512i SHUFFLE_REMOVE16 [];
#endif

#define LOCK_MASK (1ULL << 63)
#define UNLOCK_MASK ~(1ULL << 63)

static inline void lock(vqf_block& block)
{
#ifdef ENABLE_THREADS
   uint64_t *data;
#if TAG_BITS == 8
   data = block.md + 1;
#elif TAG_BITS == 16
   data = &block.md;
#endif
   while ((__sync_fetch_and_or(data, LOCK_MASK) & (1ULL << 63)) != 0) {}
#endif
}

static inline void unlock(vqf_block& block)
{
#ifdef ENABLE_THREADS
   uint64_t *data;
#if TAG_BITS == 8
   data = block.md + 1;
#elif TAG_BITS == 16
   data = &block.md;
#endif
   __sync_fetch_and_and(data, UNLOCK_MASK);
#endif
}

static inline void lock_blocks(vqf_filter * restrict filter, uint64_t index1, uint64_t index2)  {
#ifdef ENABLE_THREADS
   if (index1 < index2) {
      lock(filter->blocks[index1/QUQU_BUCKETS_PER_BLOCK]);
      lock(filter->blocks[index2/QUQU_BUCKETS_PER_BLOCK]);
   } else {
      lock(filter->blocks[index2/QUQU_BUCKETS_PER_BLOCK]);
      lock(filter->blocks[index1/QUQU_BUCKETS_PER_BLOCK]);
   }
#endif
}

static inline void unlock_blocks(vqf_filter * restrict filter, uint64_t index1, uint64_t index2)  {
#ifdef ENABLE_THREADS
   if (index1 < index2) {
      unlock(filter->blocks[index1/QUQU_BUCKETS_PER_BLOCK]);
      unlock(filter->blocks[index2/QUQU_BUCKETS_PER_BLOCK]);
   } else {
      unlock(filter->blocks[index2/QUQU_BUCKETS_PER_BLOCK]);
      unlock(filter->blocks[index1/QUQU_BUCKETS_PER_BLOCK]);
   }
#endif
}

static inline int word_rank(uint64_t val) {
   return __builtin_popcountll(val);
}

// Returns the position of the rank'th 1.  (rank = 0 returns the 1st 1)
// Returns 64 if there are fewer than rank+1 1s.
static inline uint64_t word_select(uint64_t val, int rank) {
   val = _pdep_u64(one[rank], val);
   return _tzcnt_u64(val);
}

// select(vec, 0) -> -1
// select(vec, i) -> 128, if i > popcnt(vec)
static inline int64_t select_128_old(__uint128_t vector, uint64_t rank) {
   uint64_t lower_word = vector & 0xffffffffffffffff;
   uint64_t lower_pdep = _pdep_u64(one[rank], lower_word);
   //uint64_t lower_select = word_select(lower_word, rank);
   if (lower_pdep != 0) {
      //assert(rank < word_rank(lower_word));
      return _tzcnt_u64(lower_pdep);
   }
   rank = rank - word_rank(lower_word);
   uint64_t higher_word = vector >> 64;
   return word_select(higher_word, rank) + 64;
}

static inline uint64_t lookup_64(uint64_t vector, uint64_t rank) {
   uint64_t lower_return = _pdep_u64(one[rank], vector) >> rank << (sizeof(uint64_t)/2);
   return lower_return;
}

static inline uint64_t lookup_128(uint64_t *vector, uint64_t rank) {
   uint64_t lower_word = vector[0];
   uint64_t lower_rank = word_rank(lower_word);
   uint64_t lower_return = _pdep_u64(one[rank], lower_word) >> rank << sizeof(__uint128_t);
   int64_t higher_rank = (int64_t)rank - lower_rank;
   uint64_t higher_word = vector[1];
   uint64_t higher_return = _pdep_u64(one[higher_rank], higher_word);
   higher_return <<= (64 + sizeof(__uint128_t) - rank);
   return lower_return + higher_return;
}

static inline int64_t select_64(uint64_t vector, uint64_t rank) {
   return _tzcnt_u64(lookup_64(vector, rank));
}

static inline int64_t select_128(uint64_t *vector, uint64_t rank) {
   return _tzcnt_u64(lookup_128(vector, rank));
}

//assumes little endian
#if TAG_BITS == 8
void print_bits(__uint128_t num, int numbits)
{
   int i;
   for (i = 0 ; i < numbits; i++) {
      if (i != 0 && i % 8 == 0) {
         printf(":");
      }
      printf("%d", ((num >> i) & 1) == 1);
   }
   puts("");
}

void print_tags(uint8_t *tags, uint32_t size) {
   for (uint8_t i = 0; i < size; i++)
      printf("%d ", (uint32_t)tags[i]);
   printf("\n");
}

void print_block(vqf_filter *filter, uint64_t block_index) {
   printf("block index: %ld\n", block_index);
   printf("metadata: ");
   uint64_t *md = filter->blocks[block_index].md;
   print_bits(*(__uint128_t *)md, QUQU_BUCKETS_PER_BLOCK +
         QUQU_SLOTS_PER_BLOCK);
   printf("tags: ");
   print_tags(filter->blocks[block_index].tags, QUQU_SLOTS_PER_BLOCK);
}
#elif TAG_BITS == 16
void print_bits(uint64_t num, int numbits)
{
   int i;
   for (i = 0 ; i < numbits; i++) {
      if (i != 0 && i % 8 == 0) {
         printf(":");
      }
      printf("%d", ((num >> i) & 1) == 1);
   }
   puts("");
}
void print_tags(uint16_t *tags, uint32_t size) {
   for (uint8_t i = 0; i < size; i++)
      printf("%d ", (uint32_t)tags[i]);
   printf("\n");
}
void print_block(vqf_filter *filter, uint64_t block_index) {
   printf("block index: %ld\n", block_index);
   printf("metadata: ");
   uint64_t md = filter->blocks[block_index].md;
   print_bits(md, QUQU_BUCKETS_PER_BLOCK + QUQU_SLOTS_PER_BLOCK);
   printf("tags: ");
   print_tags(filter->blocks[block_index].tags, QUQU_SLOTS_PER_BLOCK);
}
#endif

#ifdef __AVX512BW__
#if TAG_BITS == 8
static inline void update_tags_512(vqf_block * restrict block, uint8_t index, uint8_t tag) {
   block->tags[47] = tag;	// add tag at the end

   __m512i vector = _mm512_loadu_si512(reinterpret_cast<__m512i*>(block));
   vector = _mm512_permutexvar_epi8(SHUFFLE[index], vector);
   _mm512_storeu_si512(reinterpret_cast<__m512i*>(block), vector);
}

static inline void remove_tags_512(vqf_block * restrict block, uint8_t index) {
   __m512i vector = _mm512_loadu_si512(reinterpret_cast<__m512i*>(block));
   vector = _mm512_permutexvar_epi8(SHUFFLE_REMOVE[index], vector);
   _mm512_storeu_si512(reinterpret_cast<__m512i*>(block), vector);
}
#elif TAG_BITS == 16
static inline void update_tags_512(vqf_block * restrict block, uint8_t index, uint16_t tag) {
   block->tags[27] = tag;	// add tag at the end

   __m512i vector = _mm512_loadu_si512(reinterpret_cast<__m512i*>(block));
   vector = _mm512_permutexvar_epi16(SHUFFLE16[index], vector);
   _mm512_storeu_si512(reinterpret_cast<__m512i*>(block), vector);
}

static inline void remove_tags_512(vqf_block * restrict block, uint8_t index) {
   __m512i vector = _mm512_loadu_si512(reinterpret_cast<__m512i*>(block));
   vector = _mm512_permutexvar_epi16(SHUFFLE_REMOVE16[index], vector);
   _mm512_storeu_si512(reinterpret_cast<__m512i*>(block), vector);
}
#endif
#else
#if TAG_BITS == 8
static inline void update_tags_512(vqf_block * restrict block, uint8_t index, uint8_t tag) {
   index -= 16;
   memmove(&block->tags[index + 1], &block->tags[index], sizeof(block->tags) / sizeof(block->tags[0]) - index - 1);
   block->tags[index] = tag;
}

static inline void remove_tags_512(vqf_block * restrict block, uint8_t index) {
   index -= 16;
   memmove(&block->tags[index], &block->tags[index+1], sizeof(block->tags) / sizeof(block->tags[0]) - index);
}
#elif TAG_BITS == 16
static inline void update_tags_512(vqf_block * restrict block, uint8_t index, uint16_t tag) {
   index -= 4;
   memmove(&block->tags[index + 1], &block->tags[index], (sizeof(block->tags) / sizeof(block->tags[0]) - index - 1) * 2);
   block->tags[index] = tag;
}

static inline void remove_tags_512(vqf_block * restrict block, uint8_t index) {
   index -= 4;
   memmove(&block->tags[index], &block->tags[index+1], (sizeof(block->tags) / sizeof(block->tags[0]) - index) * 2);
}
#endif
#endif

#if 0
// Shuffle using AVX2 vector instruction. It turns out memmove is faster compared to AVX2.
inline __m256i cross_lane_shuffle(const __m256i & value, const __m256i &
      shuffle) 
{ 
   return _mm256_or_si256(_mm256_shuffle_epi8(value, _mm256_add_epi8(shuffle,
               K[0])), 
         _mm256_shuffle_epi8(_mm256_permute4x64_epi64(value, 0x4E),
            _mm256_add_epi8(shuffle, K[1]))); 
} 

#define SHUFFLE_SIZE 32
void shuffle_256(uint8_t * restrict source, __m256i shuffle) {
   __m256i vector = _mm256_loadu_si256(reinterpret_cast<__m256i*>(source)); 

   vector = cross_lane_shuffle(vector, shuffle); 
   _mm256_storeu_si256(reinterpret_cast<__m256i*>(source), vector); 
} 

static inline void update_tags_256(uint8_t * restrict block, uint8_t index,
      uint8_t tag) {
   index = index + sizeof(__uint128_t);	// offset index based on md field.
   block[63] = tag;	// add tag at the end
   shuffle_256(block + SHUFFLE_SIZE, RM[index]); // right block shuffle
   if (index < SHUFFLE_SIZE) {		// if index lies in the left block
      std::swap(block[31], block[32]);	// move tag to the end of left block
      shuffle_256(block, LM[index]);	// shuffle left block
   }
}
#endif

#if TAG_BITS == 8
static inline void update_md(uint64_t *md, uint8_t index) {
   uint64_t carry = (md[0] >> 63) & carry_pdep_table[index];
   md[1] = _pdep_u64(md[1],         high_order_pdep_table[index]) | carry;
   md[0] = _pdep_u64(md[0],         low_order_pdep_table[index]);
}

static inline void remove_md(uint64_t *md, uint8_t index) {
   uint64_t carry = (md[1] & carry_pdep_table[index]) << 63;
   md[1] = _pext_u64(md[1],  high_order_pdep_table[index]) | (1ULL << 63);
   md[0] = _pext_u64(md[0],  low_order_pdep_table[index]) | carry;
}

// number of 0s in the metadata is the number of tags.
static inline uint64_t get_block_free_space(uint64_t *vector) {
   uint64_t lower_word = vector[0];
   uint64_t higher_word = vector[1];
   return word_rank(lower_word) + word_rank(higher_word);
}
#elif TAG_BITS == 16
static inline void update_md(uint64_t *md, uint8_t index) {
   *md = _pdep_u64(*md, low_order_pdep_table[index]);
}

static inline void remove_md(uint64_t *md, uint8_t index) {
   *md = _pext_u64(*md, low_order_pdep_table[index]) | (1ULL << 63);
}

// number of 0s in the metadata is the number of tags.
static inline uint64_t get_block_free_space(uint64_t vector) {
   return word_rank(vector);
}
#endif

// Create n/log(n) blocks of log(n) slots.
// log(n) is 51 given a cache line size.
// n/51 blocks.
vqf_filter * vqf_init(uint64_t nslots) {
   vqf_filter *filter;

   uint64_t total_blocks = (nslots + QUQU_SLOTS_PER_BLOCK)/QUQU_SLOTS_PER_BLOCK;
   uint64_t total_size_in_bytes = sizeof(vqf_block) * total_blocks;

   filter = (vqf_filter *)malloc(sizeof(*filter) + total_size_in_bytes);
   printf("Size: %ld\n",total_size_in_bytes);
   assert(filter);

   filter->metadata.total_size_in_bytes = total_size_in_bytes;
   filter->metadata.nslots = total_blocks * QUQU_SLOTS_PER_BLOCK;
#if TAG_BITS == 8
   filter->metadata.key_remainder_bits = 8;
#elif TAG_BITS == 16
   filter->metadata.key_remainder_bits = 16;
#endif
   filter->metadata.range = total_blocks * QUQU_BUCKETS_PER_BLOCK * (1ULL << filter->metadata.key_remainder_bits);
   filter->metadata.nblocks = total_blocks;
   filter->metadata.nelts = 0;

   // memset to 1
#if TAG_BITS == 8
   for (uint64_t i = 0; i < total_blocks; i++) {
      filter->blocks[i].md[0] = UINT64_MAX;
      filter->blocks[i].md[1] = UINT64_MAX;
      // reset the most significant bit of metadata for locking.
      filter->blocks[i].md[1] = filter->blocks[i].md[1] & ~(1ULL << 63);
   }
#elif TAG_BITS == 16
   for (uint64_t i = 0; i < total_blocks; i++) {
      filter->blocks[i].md = UINT64_MAX;
      filter->blocks[i].md = filter->blocks[i].md & ~(1ULL << 63);
   }
#endif

   return filter;
}

// If the item goes in the i'th slot (starting from 0) in the block then
// find the i'th 0 in the metadata, insert a 1 after that and shift the rest
// by 1 bit.
// Insert the new tag at the end of its run and shift the rest by 1 slot.
bool vqf_insert(vqf_filter * restrict filter, uint64_t hash) {
   vqf_metadata * restrict metadata           = &filter->metadata;
   vqf_block    * restrict blocks             = filter->blocks;
   uint64_t                 key_remainder_bits = metadata->key_remainder_bits;
   uint64_t                 range              = metadata->range;

   uint64_t block_index = hash >> key_remainder_bits;
   lock(blocks[block_index/QUQU_BUCKETS_PER_BLOCK]);
#if TAG_BITS == 8
   uint64_t *block_md = blocks[block_index/QUQU_BUCKETS_PER_BLOCK].md;
   uint64_t block_free = get_block_free_space(block_md);
#elif TAG_BITS == 16
   uint64_t *block_md = &blocks[block_index/QUQU_BUCKETS_PER_BLOCK].md;
   uint64_t block_free = get_block_free_space(*block_md);
#endif
   uint64_t tag = hash & TAG_MASK;
   uint64_t alt_block_index = ((hash ^ (tag * 0x5bd1e995)) % range) >> key_remainder_bits;

   __builtin_prefetch(&blocks[alt_block_index/QUQU_BUCKETS_PER_BLOCK]);

   if (block_free < QUQU_CHECK_ALT && block_index/QUQU_BUCKETS_PER_BLOCK != alt_block_index/QUQU_BUCKETS_PER_BLOCK) {
      unlock(blocks[block_index/QUQU_BUCKETS_PER_BLOCK]);
      lock_blocks(filter, block_index, alt_block_index);
#if TAG_BITS == 8
      uint64_t *alt_block_md = blocks[alt_block_index/QUQU_BUCKETS_PER_BLOCK].md;
      uint64_t alt_block_free = get_block_free_space(alt_block_md);
#elif TAG_BITS == 16
      uint64_t *alt_block_md = &blocks[alt_block_index/QUQU_BUCKETS_PER_BLOCK].md;
      uint64_t alt_block_free = get_block_free_space(*alt_block_md);
#endif
      // pick the least loaded block
      if (alt_block_free > block_free) {
         unlock(blocks[block_index/QUQU_BUCKETS_PER_BLOCK]);
         block_index = alt_block_index;
         block_md = alt_block_md;
      } else if (block_free == QUQU_BUCKETS_PER_BLOCK) {
         unlock_blocks(filter, block_index, alt_block_index);
         fprintf(stderr, "vqf filter is full.");
         return false;
         //exit(EXIT_FAILURE);
      } else {
         unlock(blocks[alt_block_index/QUQU_BUCKETS_PER_BLOCK]);
      }

   }

   uint64_t index = block_index / QUQU_BUCKETS_PER_BLOCK;
   uint64_t offset = block_index % QUQU_BUCKETS_PER_BLOCK;

#if TAG_BITS == 8
   uint64_t slot_index = select_128(block_md, offset);
   uint64_t select_index = slot_index + offset - sizeof(__uint128_t);
#elif TAG_BITS == 16
   uint64_t slot_index = select_64(*block_md, offset);
   uint64_t select_index = slot_index + offset - (sizeof(uint64_t)/2);
#endif
   /*printf("index: %ld tag: %ld offset: %ld\n", index, tag, offset);*/
   /*print_block(filter, index);*/

   update_tags_512(&blocks[index], slot_index,tag);
   update_md(block_md, select_index);
   /*print_block(filter, index);*/
   unlock(blocks[block_index/QUQU_BUCKETS_PER_BLOCK]);
   return true;
}

static inline bool remove_tags(vqf_filter * restrict filter, uint64_t tag,
      uint64_t block_index) {
   uint64_t index = block_index / QUQU_BUCKETS_PER_BLOCK;
   uint64_t offset = block_index % QUQU_BUCKETS_PER_BLOCK;

#ifdef __AVX512BW__
#if TAG_BITS == 8
   __m512i bcast = _mm512_set1_epi8(tag);
   __m512i block =
      _mm512_loadu_si512(reinterpret_cast<__m512i*>(&filter->blocks[index]));
   volatile __mmask64 result = _mm512_cmp_epi8_mask(bcast, block, _MM_CMPINT_EQ);
#elif TAG_BITS == 16
   __m512i bcast = _mm512_set1_epi16(tag);
   __m512i block =
      _mm512_loadu_si512(reinterpret_cast<__m512i*>(&filter->blocks[index]));
   volatile __mmask64 result = _mm512_cmp_epi16_mask(bcast, block, _MM_CMPINT_EQ);
#endif
#else
#if TAG_BITS == 8
   __m256i bcast = _mm256_set1_epi8(tag);
   __m256i block = _mm256_loadu_si256(reinterpret_cast<__m256i*>(&filter->blocks[index]));
   __m256i result1t = _mm256_cmpeq_epi8(bcast, block);
   __mmask32 result1 = _mm256_movemask_epi8(result1t);
   /*__mmask32 result1 = _mm256_cmp_epi8_mask(bcast, block, _MM_CMPINT_EQ);*/
   block = _mm256_loadu_si256(reinterpret_cast<__m256i*>((uint8_t*)&filter->blocks[index]+32));
   __m256i result2t = _mm256_cmpeq_epi8(bcast, block);
   __mmask32 result2 = _mm256_movemask_epi8(result2t);
   /*__mmask32 result2 = _mm256_cmp_epi8_mask(bcast, block, _MM_CMPINT_EQ);*/
   uint64_t result = (uint64_t)result2 << 32 | (uint64_t)result1;
#elif TAG_BITS == 16
   uint64_t alt_mask = 0x55555555;
   __m256i bcast = _mm256_set1_epi16(tag);
   __m256i block = _mm256_loadu_si256(reinterpret_cast<__m256i*>(&filter->blocks[index]));
   __m256i result1t = _mm256_cmpeq_epi16(bcast, block);
   __mmask32 result1 = _mm256_movemask_epi8(result1t);
   result1 = _pext_u32(result1, alt_mask);
   /*__mmask32 result1 = _mm256_cmp_epi8_mask(bcast, block, _MM_CMPINT_EQ);*/
   block = _mm256_loadu_si256(reinterpret_cast<__m256i*>((uint8_t*)&filter->blocks[index]+32));
   __m256i result2t = _mm256_cmpeq_epi16(bcast, block);
   __mmask32 result2 = _mm256_movemask_epi8(result2t);
   result2 = _pext_u32(result2, alt_mask);
   /*__mmask32 result2 = _mm256_cmp_epi8_mask(bcast, block, _MM_CMPINT_EQ);*/
   uint64_t result = (uint64_t)result2 << 16 | (uint64_t)result1;
#endif
#endif

   if (result == 0) {
      // no matching tags, can bail
      return false;
   }

#if TAG_BITS == 8
   uint64_t start = offset != 0 ? lookup_128(filter->blocks[index].md, offset -
         1) : one[0] << 2 * sizeof(uint64_t);
   uint64_t end = lookup_128(filter->blocks[index].md, offset);
#elif TAG_BITS == 16
   uint64_t start = offset != 0 ? lookup_64(filter->blocks[index].md, offset -
         1) : one[0] << (sizeof(uint64_t)/2);
   uint64_t end = lookup_64(filter->blocks[index].md, offset);
#endif
   uint64_t mask = end - start;

   uint64_t check_indexes = mask & result;
   if (check_indexes != 0) { // remove the first available tag
      vqf_block    * restrict blocks             = filter->blocks;
      uint64_t remove_index = __builtin_ctzll(check_indexes);
      remove_tags_512(&blocks[index], remove_index);
#if TAG_BITS == 8
      remove_index = remove_index + offset - sizeof(__uint128_t);
      uint64_t *block_md = blocks[block_index / QUQU_BUCKETS_PER_BLOCK].md;
      remove_md(block_md, remove_index);
#elif TAG_BITS == 16
      remove_index = remove_index + offset - sizeof(uint64_t);
      uint64_t *block_md = &blocks[block_index / QUQU_BUCKETS_PER_BLOCK].md;
      remove_md(block_md, remove_index);
#endif
      return true;
   } else
      return false;
}

bool vqf_remove(vqf_filter * restrict filter, uint64_t hash) {
   vqf_metadata * restrict metadata           = &filter->metadata;
   uint64_t                 key_remainder_bits = metadata->key_remainder_bits;
   uint64_t                 range              = metadata->range;

   uint64_t block_index = hash >> key_remainder_bits;
   uint64_t tag = hash & TAG_MASK;
   uint64_t alt_block_index = ((hash ^ (tag * 0x5bd1e995)) % range) >> key_remainder_bits;

   __builtin_prefetch(&filter->blocks[alt_block_index / QUQU_BUCKETS_PER_BLOCK]);

   return remove_tags(filter, tag, block_index) || remove_tags(filter, tag, alt_block_index);
}

static inline bool check_tags(vqf_filter * restrict filter, uint64_t tag,
      uint64_t block_index) {
   uint64_t index = block_index / QUQU_BUCKETS_PER_BLOCK;
   uint64_t offset = block_index % QUQU_BUCKETS_PER_BLOCK;

#ifdef __AVX512BW__
#if TAG_BITS == 8
   __m512i bcast = _mm512_set1_epi8(tag);
   __m512i block =
      _mm512_loadu_si512(reinterpret_cast<__m512i*>(&filter->blocks[index]));
   volatile __mmask64 result = _mm512_cmp_epi8_mask(bcast, block, _MM_CMPINT_EQ);
#elif TAG_BITS == 16
   __m512i bcast = _mm512_set1_epi16(tag);
   __m512i block =
      _mm512_loadu_si512(reinterpret_cast<__m512i*>(&filter->blocks[index]));
   volatile __mmask64 result = _mm512_cmp_epi16_mask(bcast, block, _MM_CMPINT_EQ);
#endif
#else
#if TAG_BITS == 8
   __m256i bcast = _mm256_set1_epi8(tag);
   __m256i block = _mm256_loadu_si256(reinterpret_cast<__m256i*>(&filter->blocks[index]));
   __m256i result1t = _mm256_cmpeq_epi8(bcast, block);
   __mmask32 result1 = _mm256_movemask_epi8(result1t);
   /*__mmask32 result1 = _mm256_cmp_epi8_mask(bcast, block, _MM_CMPINT_EQ);*/
   block = _mm256_loadu_si256(reinterpret_cast<__m256i*>((uint8_t*)&filter->blocks[index]+32));
   __m256i result2t = _mm256_cmpeq_epi8(bcast, block);
   __mmask32 result2 = _mm256_movemask_epi8(result2t);
   /*__mmask32 result2 = _mm256_cmp_epi8_mask(bcast, block, _MM_CMPINT_EQ);*/
   uint64_t result = (uint64_t)result2 << 32 | (uint64_t)result1;
#elif TAG_BITS == 16
   uint64_t alt_mask = 0x55555555;
   __m256i bcast = _mm256_set1_epi16(tag);
   __m256i block = _mm256_loadu_si256(reinterpret_cast<__m256i*>(&filter->blocks[index]));
   __m256i result1t = _mm256_cmpeq_epi16(bcast, block);
   __mmask32 result1 = _mm256_movemask_epi8(result1t);
   result1 = _pext_u32(result1, alt_mask);
   /*__mmask32 result1 = _mm256_cmp_epi8_mask(bcast, block, _MM_CMPINT_EQ);*/
   block = _mm256_loadu_si256(reinterpret_cast<__m256i*>((uint8_t*)&filter->blocks[index]+32));
   __m256i result2t = _mm256_cmpeq_epi16(bcast, block);
   __mmask32 result2 = _mm256_movemask_epi8(result2t);
   result2 = _pext_u32(result2, alt_mask);
   /*__mmask32 result2 = _mm256_cmp_epi8_mask(bcast, block, _MM_CMPINT_EQ);*/
   uint64_t result = (uint64_t)result2 << 16 | (uint64_t)result1;
#endif
#endif

   if (result == 0) {
      // no matching tags, can bail
      return false;
   }

#if TAG_BITS == 8
   uint64_t start = offset != 0 ? lookup_128(filter->blocks[index].md, offset -
         1) : one[0] << 2 * sizeof(uint64_t);
   uint64_t end = lookup_128(filter->blocks[index].md, offset);
#elif TAG_BITS == 16
   uint64_t start = offset != 0 ? lookup_64(filter->blocks[index].md, offset -
         1) : one[0] << (sizeof(uint64_t)/2);
   uint64_t end = lookup_64(filter->blocks[index].md, offset);
#endif
   uint64_t mask = end - start;
   return (mask & result) != 0;
}

// If the item goes in the i'th slot (starting from 0) in the block then
// select(i) - i is the slot index for the end of the run.
bool vqf_is_present(vqf_filter * restrict filter, uint64_t hash) {
   vqf_metadata * restrict metadata           = &filter->metadata;
   //vqf_block    * restrict blocks             = filter->blocks;
   uint64_t                 key_remainder_bits = metadata->key_remainder_bits;
   uint64_t                 range              = metadata->range;

   uint64_t block_index = hash >> key_remainder_bits;
   uint64_t tag = hash & TAG_MASK;
   uint64_t alt_block_index = ((hash ^ (tag * 0x5bd1e995)) % range) >> key_remainder_bits;

   __builtin_prefetch(&filter->blocks[alt_block_index / QUQU_BUCKETS_PER_BLOCK]);

   return check_tags(filter, tag, block_index) || check_tags(filter, tag, alt_block_index);

   /*if (!ret) {*/
   /*printf("tag: %ld offset: %ld\n", tag, block_index % QUQU_SLOTS_PER_BLOCK);*/
   /*print_block(filter, block_index / QUQU_SLOTS_PER_BLOCK);*/
   /*print_block(filter, alt_block_index / QUQU_SLOTS_PER_BLOCK);*/
   /*}*/
}

} // namespace

#endif
