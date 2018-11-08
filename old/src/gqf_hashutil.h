/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *
 * ============================================================================
 */

#ifndef _HASHUTIL_H_
#define _HASHUTIL_H_

#include <sys/types.h>
#include <stdlib.h>
#include <stdint.h>

uint64_t MurmurHash64B ( const void * key, int len, unsigned int seed );
uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed );

uint64_t hash_64(uint64_t key, uint64_t mask);
uint64_t hash_64i(uint64_t key, uint64_t mask);

#endif  // #ifndef _HASHUTIL_H_


