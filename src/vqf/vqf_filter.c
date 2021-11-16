#ifndef _VQF_FILTER_8BIT_C_
#define _VQF_FILTER_8BIT_C_

// https://github.com/splatlab/vqf

/*
 * ============================================================================
 *
 *       Filename:  vqf_filter_8bit.c
 *
 *         Author:  Prashant Pandey (), ppandey@berkeley.edu, LBNL/UCB
 *       Author 2:  Thomas Mueller Graf
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

#define TAG_BITS 8
#define TAG_MASK 0xff
#define QUQU_SLOTS_PER_BLOCK 48
#define QUQU_BUCKETS_PER_BLOCK 80
#define QUQU_CHECK_ALT 92

static inline __attribute__((always_inline)) uint64_t rotateLeft(uint64_t n, unsigned int c) {
    const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
    c &= mask;
    return (n << c) | (n >> ((-c) & mask));
}


static inline int word_rank(uint64_t val) {
   return __builtin_popcountll(val);
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

static inline int64_t select_128(uint64_t *vector, uint64_t rank) {
   return _tzcnt_u64(lookup_128(vector, rank));
}

static inline void update_tags_512(vqf_block * restrict block, uint8_t index, uint8_t tag) {
   index -= 16;
   memmove(&block->tags[index + 1], &block->tags[index], sizeof(block->tags) / sizeof(block->tags[0]) - index - 1);
   block->tags[index] = tag;
}

static inline void update_md(uint64_t *md, uint8_t index) {
   uint64_t carry = (md[0] >> 63) & carry_pdep_table[index];
   md[1] = _pdep_u64(md[1],         high_order_pdep_table[index]) | carry;
   md[0] = _pdep_u64(md[0],         low_order_pdep_table[index]);
}

// number of 0s in the metadata is the number of tags.
static inline uint64_t get_block_free_space(uint64_t *vector) {
   uint64_t lower_word = vector[0];
   uint64_t higher_word = vector[1];
   return word_rank(lower_word) + word_rank(higher_word);
}

// Create n/log(n) blocks of log(n) slots.
// log(n) is 51 given a cache line size.
// n/51 blocks.
vqf_filter * vqf_init(uint64_t nslots) {
   vqf_filter *filter;

   uint64_t total_blocks = (nslots + QUQU_SLOTS_PER_BLOCK)/QUQU_SLOTS_PER_BLOCK;
   uint64_t total_size_in_bytes = sizeof(vqf_block) * total_blocks;

   filter = (vqf_filter *)malloc(sizeof(*filter) + total_size_in_bytes);
   // printf("Size: %ld\n",total_size_in_bytes);
   assert(filter);

   filter->metadata.total_size_in_bytes = total_size_in_bytes;
   filter->metadata.nslots = total_blocks * QUQU_SLOTS_PER_BLOCK;
   filter->metadata.range = total_blocks * QUQU_BUCKETS_PER_BLOCK;
   filter->metadata.nblocks = total_blocks;
   filter->metadata.nelts = 0;

   // memset to 1
   for (uint64_t i = 0; i < total_blocks; i++) {
      filter->blocks[i].md[0] = UINT64_MAX;
      filter->blocks[i].md[1] = UINT64_MAX;
      // reset the most significant bit of metadata for locking.
      filter->blocks[i].md[1] = filter->blocks[i].md[1] & ~(1ULL << 63);
   }

   return filter;
}

static inline bool check_tags(vqf_filter * restrict filter, uint64_t tag,
      uint64_t block_index) {
   uint64_t index = block_index / QUQU_BUCKETS_PER_BLOCK;
   uint64_t offset = block_index % QUQU_BUCKETS_PER_BLOCK;

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

   if (result == 0) {
      // no matching tags, can bail
      return false;
   }

   uint64_t start = offset != 0 ? lookup_128(filter->blocks[index].md, offset -
         1) : one[0] << 2 * sizeof(uint64_t);
   uint64_t end = lookup_128(filter->blocks[index].md, offset);
   uint64_t mask = end - start;
   return (mask & result) != 0;
}

// If the item goes in the i'th slot (starting from 0) in the block then
// select(i) - i is the slot index for the end of the run.
bool vqf_is_present(vqf_filter * restrict filter, uint64_t h) {

    vqf_metadata * restrict metadata           = &filter->metadata;
    //vqf_block    * restrict blocks             = filter->blocks;
    uint64_t                 range              = metadata->range;
    uint64_t h2 = rotateLeft(h, 32);
    __uint128_t x = (__uint128_t)h * (__uint128_t)range;
    uint64_t block_index = (uint64_t)(x >> 64);
    __uint128_t y = (__uint128_t)h2 * (__uint128_t)range;
    uint64_t alt_block_index = (uint64_t)(y >> 64);
    uint64_t tag = (h ^ h2) & TAG_MASK;

    __builtin_prefetch(&filter->blocks[alt_block_index / QUQU_BUCKETS_PER_BLOCK]);

   return check_tags(filter, tag, block_index) || check_tags(filter, tag, alt_block_index);

}

// If the item goes in the i'th slot (starting from 0) in the block then
// find the i'th 0 in the metadata, insert a 1 after that and shift the rest
// by 1 bit.
// Insert the new tag at the end of its run and shift the rest by 1 slot.
bool vqf_insert(vqf_filter * restrict filter, uint64_t h) {

   vqf_metadata * restrict metadata           = &filter->metadata;
   vqf_block    * restrict blocks             = filter->blocks;
   uint64_t                 range              = metadata->range;
    uint64_t h2 = rotateLeft(h, 32);
    __uint128_t x = (__uint128_t)h * (__uint128_t)range;
    uint64_t block_index = (uint64_t)(x >> 64);
    __uint128_t y = (__uint128_t)h2 * (__uint128_t)range;
    uint64_t alt_block_index = (uint64_t)(y >> 64);
    uint64_t tag = (h ^ h2) & TAG_MASK;

    uint64_t *block_md = blocks[block_index/QUQU_BUCKETS_PER_BLOCK].md;
    uint64_t block_free = get_block_free_space(block_md);

    __builtin_prefetch(&blocks[alt_block_index/QUQU_BUCKETS_PER_BLOCK]);

    if (block_free < QUQU_CHECK_ALT && block_index/QUQU_BUCKETS_PER_BLOCK != alt_block_index/QUQU_BUCKETS_PER_BLOCK) {
        uint64_t *alt_block_md = blocks[alt_block_index/QUQU_BUCKETS_PER_BLOCK].md;
        uint64_t alt_block_free = get_block_free_space(alt_block_md);
        // pick the least loaded block
        if (alt_block_free > block_free) {
            block_index = alt_block_index;
            block_md = alt_block_md;
        } else if (block_free == QUQU_BUCKETS_PER_BLOCK) {
            fprintf(stderr, "vqf filter is full.");
            return false;
            //exit(EXIT_FAILURE);
        }
    }

    uint64_t index = block_index / QUQU_BUCKETS_PER_BLOCK;
    uint64_t offset = block_index % QUQU_BUCKETS_PER_BLOCK;

    uint64_t slot_index = select_128(block_md, offset);
    uint64_t select_index = slot_index + offset - sizeof(__uint128_t);

    update_tags_512(&blocks[index], slot_index,tag);
    update_md(block_md, select_index);
    return true;
}

} // namespace

#endif
