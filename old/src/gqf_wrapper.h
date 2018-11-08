/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *
 * ============================================================================
 */

#ifndef GQF_WRAPPER_H
#define GQF_WRAPPER_H

#include "gqf.h"
#include "gqf_int.h"
#include "gqf_file.h"

QF g_quotient_filter;
QFi g_quotient_filter_itr;

extern inline int gqf_init(uint64_t nbits, uint64_t num_hash_bits)
{
	uint64_t nslots = 1 << nbits;
	qf_malloc(&g_quotient_filter, nslots, num_hash_bits, 0, QF_HASH_NONE, 0);
	return 0;
}

extern inline int gqf_insert(__uint128_t val, uint64_t count)
{
	qf_insert(&g_quotient_filter, val, 0, count, QF_NO_LOCK);
	return 0;
}

extern inline int gqf_lookup(__uint128_t val)
{
	return qf_count_key_value(&g_quotient_filter, val, 0, 0);
}

extern inline __uint128_t gqf_range()
{
	return g_quotient_filter.metadata->range;
}

extern inline int gqf_destroy()
{
	qf_free(&g_quotient_filter);
	return 0;
}

extern inline int gqf_iterator(uint64_t pos)
{
	qf_iterator_from_position(&g_quotient_filter, &g_quotient_filter_itr, pos);
	return 0;
}

/* Returns 0 if the iterator is still valid (i.e. has not reached the
 * end of the QF. */
extern inline int gqf_get(uint64_t *key, uint64_t *value, uint64_t *count)
{
	return qfi_get_hash(&g_quotient_filter_itr, key, value, count);
}

/* Advance to next entry.  Returns whether or not another entry is
 * found.  */
extern inline int gqf_next()
{
	return qfi_next(&g_quotient_filter_itr);
}

/* Check to see if the if the end of the QF */
extern inline int gqf_end()
{
	return qfi_end(&g_quotient_filter_itr);
}

#endif
