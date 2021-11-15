#ifndef _VQF_FILTER_H_
#define _VQF_FILTER_H_

// https://github.com/splatlab/vqf

/*
 * ============================================================================
 *
 *       Filename:  vqf_filter.h
 *
 *         Author:  Prashant Pandey (), ppandey@berkeley.edu
 *   Organization: 	LBNL/UCB 
 *
 * ============================================================================
 */

#include <inttypes.h>
#include <stdbool.h>

#ifdef __cplusplus
#define restrict __restrict__
extern "C" {
#endif

namespace vqfilter {

	// metadata: 1 --> end of the run
	// Each 1 is preceded by k 0s, where k is the number of remainders in that
	// run.

	// We are using 8-bit tags.
	// One block consists of 48 8-bit slots covering 80 buckets, and 80+48 = 128
	// bits of metadata.
	typedef struct __attribute__ ((__packed__)) vqf_block {
		uint64_t md[2];
		uint8_t tags[48];
	} vqf_block;

	typedef struct vqf_metadata {
		uint64_t total_size_in_bytes;
		uint64_t key_remainder_bits;
		uint64_t range;
		uint64_t nblocks;
		uint64_t nelts;
		uint64_t nslots;
	} vqf_metadata;

	typedef struct vqf_filter {
		vqf_metadata metadata;
		vqf_block blocks[];
	} vqf_filter;

	vqf_filter * vqf_init(uint64_t nslots);

	bool vqf_insert(vqf_filter * restrict filter, uint64_t hash);
	
	bool vqf_remove(vqf_filter * restrict filter, uint64_t hash);

	bool vqf_is_present(vqf_filter * restrict filter, uint64_t hash);

}

#ifdef __cplusplus
}
#endif

#endif // _VQF_FILTER_H_


