/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>
 *
 * ============================================================================
 */

#include <stdlib.h>
# include <assert.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

// #include "gqf_hashutil.h"
// #include "gqf.h"
// #include "gqf_int.h"
// #include "gqf.c"

/******************************************************************
 * Code for managing the metadata bits and slots w/o interpreting *
 * the content of the slots.
 ******************************************************************/

#define MAX_VALUE(nbits) ((1ULL << (nbits)) - 1)
#define BITMASK(nbits)                                    \
  ((nbits) == 64 ? 0xffffffffffffffff : MAX_VALUE(nbits))
#define NUM_SLOTS_TO_LOCK (1ULL<<16)
#define CLUSTER_SIZE (1ULL<<14)
#define METADATA_WORD(qf,field,slot_index)                              \
  (get_block((qf), (slot_index) /                                       \
             QF_SLOTS_PER_BLOCK)->field[((slot_index)  % QF_SLOTS_PER_BLOCK) / 64])

#define GET_NO_LOCK(flag) (flag & QF_NO_LOCK)
#define GET_TRY_ONCE_LOCK(flag) (flag & QF_TRY_ONCE_LOCK)
#define GET_WAIT_FOR_LOCK(flag) (flag & QF_WAIT_FOR_LOCK)
#define GET_KEY_HASH(flag) (flag & QF_KEY_IS_HASH)

#define DISTANCE_FROM_HOME_SLOT_CUTOFF 1000
#define BILLION 1000000000L

#ifdef DEBUG
#define PRINT_DEBUG 1
#else
#define PRINT_DEBUG 0
#endif

#define DEBUG_CQF(fmt, ...) \
	do { if (PRINT_DEBUG) fprintf(stderr, fmt, __VA_ARGS__); } while (0)

#define DEBUG_DUMP(qf) \
	do { if (PRINT_DEBUG) qf_dump_metadata(qf); } while (0)

// static __inline__ unsigned long long rdtsc(void)
// {
// 	unsigned hi, lo;
// 	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
// 	return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
// }

#ifdef LOG_WAIT_TIME
static inline bool qf_spin_lock(QF *qf, volatile int *lock, uint64_t idx,
																uint8_t flag)
{
	struct timespec start, end;
	bool ret;

	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
	if (GET_WAIT_FOR_LOCK(flag) != QF_WAIT_FOR_LOCK) {
		ret = !__sync_lock_test_and_set(lock, 1);
		clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
		qf->runtimedata->wait_times[idx].locks_acquired_single_attempt++;
		qf->runtimedata->wait_times[idx].total_time_single += BILLION * (end.tv_sec -
																												start.tv_sec) +
			end.tv_nsec - start.tv_nsec;
	} else {
		if (!__sync_lock_test_and_set(lock, 1)) {
			clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
			qf->runtimedata->wait_times[idx].locks_acquired_single_attempt++;
			qf->runtimedata->wait_times[idx].total_time_single += BILLION * (end.tv_sec -
																													start.tv_sec) +
			end.tv_nsec - start.tv_nsec;
		} else {
			while (__sync_lock_test_and_set(lock, 1))
				while (*lock);
			clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
			qf->runtimedata->wait_times[idx].total_time_spinning += BILLION * (end.tv_sec -
																														start.tv_sec) +
				end.tv_nsec - start.tv_nsec;
		}
		ret = true;
	}
	qf->runtimedata->wait_times[idx].locks_taken++;

	return ret;

	/*start = rdtsc();*/
	/*if (!__sync_lock_test_and_set(lock, 1)) {*/
		/*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);*/
		/*qf->runtimedata->wait_times[idx].locks_acquired_single_attempt++;*/
		/*qf->runtimedata->wait_times[idx].total_time_single += BILLION * (end.tv_sec -
		 * start.tv_sec) + end.tv_nsec - start.tv_nsec;*/
	/*} else {*/
		/*while (__sync_lock_test_and_set(lock, 1))*/
			/*while (*lock);*/
		/*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);*/
		/*qf->runtimedata->wait_times[idx].total_time_spinning += BILLION * (end.tv_sec -
		 * start.tv_sec) + end.tv_nsec - start.tv_nsec;*/
	/*}*/

	/*end = rdtsc();*/
	/*qf->runtimedata->wait_times[idx].locks_taken++;*/
	/*return;*/
}
#else
/**
 * Try to acquire a lock once and return even if the lock is busy.
 * If spin flag is set, then spin until the lock is available.
 */
static inline bool qf_spin_lock(volatile int *lock, uint8_t flag)
{
	if (GET_WAIT_FOR_LOCK(flag) != QF_WAIT_FOR_LOCK) {
		return !__sync_lock_test_and_set(lock, 1);
	} else {
		while (__sync_lock_test_and_set(lock, 1))
			while (*lock);
		return true;
	}

	return false;
}
#endif

static inline void qf_spin_unlock(volatile int *lock)
{
	__sync_lock_release(lock);
	return;
}

static bool qf_lock(QF *qf, uint64_t hash_bucket_index, bool small, uint8_t
										runtime_lock)
{
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;
	if (small) {
#ifdef LOG_WAIT_TIME
		if (!qf_spin_lock(qf, &qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK],
											hash_bucket_index/NUM_SLOTS_TO_LOCK,
											runtime_lock))
			return false;
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			if (!qf_spin_lock(qf, &qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
												hash_bucket_index/NUM_SLOTS_TO_LOCK+1,
												runtime_lock)) {
				qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
#else
		if (!qf_spin_lock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK],
											runtime_lock))
			return false;
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			if (!qf_spin_lock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
												runtime_lock)) {
				qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
#endif
	} else {
#ifdef LOG_WAIT_TIME
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
				CLUSTER_SIZE) {
			if (!qf_spin_lock(qf,
												&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1],
												runtime_lock))
				return false;
		}
		if (!qf_spin_lock(qf,
											&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK],
											runtime_lock)) {
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
		if (!qf_spin_lock(qf, &qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
											runtime_lock)) {
			qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
#else
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
				CLUSTER_SIZE) {
			if
				(!qf_spin_lock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1],
											 runtime_lock))
				return false;
		}
		if (!qf_spin_lock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK],
											runtime_lock)) {
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
		if (!qf_spin_lock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
											runtime_lock)) {
			qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
#endif
	}
	return true;
}

static void qf_unlock(QF *qf, uint64_t hash_bucket_index, bool small)
{
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;
	if (small) {
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
		}
		qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
	} else {
		qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
		qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
				CLUSTER_SIZE)
			qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
	}
}

static void modify_metadata(QF *qf, uint64_t *metadata, int cnt)
{
#ifdef LOG_WAIT_TIME
	qf_spin_lock(qf, &qf->runtimedata->metadata_lock,
							 qf->runtimedata->num_locks, QF_WAIT_FOR_LOCK);
#else
	qf_spin_lock(&qf->runtimedata->metadata_lock, QF_WAIT_FOR_LOCK);
#endif
	*metadata = *metadata + cnt;
	qf_spin_unlock(&qf->runtimedata->metadata_lock);
	return;
}

static inline int popcnt(uint64_t val)
{
	asm("popcnt %[val], %[val]"
			: [val] "+r" (val)
			:
			: "cc");
	return val;
}

static inline int64_t bitscanreverse(uint64_t val)
{
	if (val == 0) {
		return -1;
	} else {
		asm("bsr %[val], %[val]"
				: [val] "+r" (val)
				:
				: "cc");
		return val;
	}
}

static inline int popcntv(const uint64_t val, int ignore)
{
	if (ignore % 64)
		return popcnt (val & ~BITMASK(ignore % 64));
	else
		return popcnt(val);
}

// Returns the number of 1s up to (and including) the pos'th bit
// Bits are numbered from 0
static inline int bitrank(uint64_t val, int pos) {
	val = val & ((2ULL << pos) - 1);
	asm("popcnt %[val], %[val]"
			: [val] "+r" (val)
			:
			: "cc");
	return val;
}

/**
 * Returns the position of the k-th 1 in the 64-bit word x.
 * k is 0-based, so k=0 returns the position of the first 1.
 *
 * Uses the broadword selection algorithm by Vigna [1], improved by Gog
 * and Petri [2] and Vigna [3].
 *
 * [1] Sebastiano Vigna. Broadword Implementation of Rank/Select
 *    Queries. WEA, 2008
 *
 * [2] Simon Gog, Matthias Petri. Optimized succinct data
 * structures for massive data. Softw. Pract. Exper., 2014
 *
 * [3] Sebastiano Vigna. MG4J 5.2.1. http://mg4j.di.unimi.it/
 * The following code is taken from
 * https://github.com/facebook/folly/blob/b28186247104f8b90cfbe094d289c91f9e413317/folly/experimental/Select64.h
 */
const uint8_t kSelectInByte[2048] = {
	8, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0,
	1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0,
	2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0,
	1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0,
	3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7, 0,
	1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0,
	2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0,
	1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0,
	1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8, 8, 8, 1,
	8, 2, 2, 1, 8, 3, 3, 1, 3, 2, 2, 1, 8, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2,
	2, 1, 8, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1,
	4, 3, 3, 1, 3, 2, 2, 1, 8, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4,
	4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1,
	3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 7, 7, 1, 7, 2,
	2, 1, 7, 3, 3, 1, 3, 2, 2, 1, 7, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
	7, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3,
	3, 1, 3, 2, 2, 1, 7, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1,
	4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2,
	2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 8, 8, 8, 8, 8, 8, 2,
	8, 8, 8, 3, 8, 3, 3, 2, 8, 8, 8, 4, 8, 4, 4, 2, 8, 4, 4, 3, 4, 3, 3, 2, 8, 8,
	8, 5, 8, 5, 5, 2, 8, 5, 5, 3, 5, 3, 3, 2, 8, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3,
	4, 3, 3, 2, 8, 8, 8, 6, 8, 6, 6, 2, 8, 6, 6, 3, 6, 3, 3, 2, 8, 6, 6, 4, 6, 4,
	4, 2, 6, 4, 4, 3, 4, 3, 3, 2, 8, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2,
	6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 7, 8, 7, 7, 2, 8, 7,
	7, 3, 7, 3, 3, 2, 8, 7, 7, 4, 7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2, 8, 7, 7, 5,
	7, 5, 5, 2, 7, 5, 5, 3, 5, 3, 3, 2, 7, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3,
	3, 2, 8, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2, 7, 6, 6, 4, 6, 4, 4, 2,
	6, 4, 4, 3, 4, 3, 3, 2, 7, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2, 6, 5,
	5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 4, 8, 4, 4, 3, 8, 8, 8, 8, 8, 8,
	8, 5, 8, 8, 8, 5, 8, 5, 5, 3, 8, 8, 8, 5, 8, 5, 5, 4, 8, 5, 5, 4, 5, 4, 4, 3,
	8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 3, 8, 8, 8, 6, 8, 6, 6, 4, 8, 6,
	6, 4, 6, 4, 4, 3, 8, 8, 8, 6, 8, 6, 6, 5, 8, 6, 6, 5, 6, 5, 5, 3, 8, 6, 6, 5,
	6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7,
	7, 3, 8, 8, 8, 7, 8, 7, 7, 4, 8, 7, 7, 4, 7, 4, 4, 3, 8, 8, 8, 7, 8, 7, 7, 5,
	8, 7, 7, 5, 7, 5, 5, 3, 8, 7, 7, 5, 7, 5, 5, 4, 7, 5, 5, 4, 5, 4, 4, 3, 8, 8,
	8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 3, 8, 7, 7, 6, 7, 6, 6, 4, 7, 6, 6, 4,
	6, 4, 4, 3, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 3, 7, 6, 6, 5, 6, 5,
	5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 5, 8, 5, 5, 4, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6,
	6, 4, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 5, 8, 8, 8, 6, 8, 6, 6, 5,
	8, 6, 6, 5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8,
	8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 4, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7,
	8, 7, 7, 5, 8, 8, 8, 7, 8, 7, 7, 5, 8, 7, 7, 5, 7, 5, 5, 4, 8, 8, 8, 8, 8, 8,
	8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 4,
	8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 5, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6,
	6, 5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6,
	8, 6, 6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7,
	8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 8,
	8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6,
	6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7
};

static inline uint64_t _select64(uint64_t x, int k)
{
	if (k >= popcnt(x)) { return 64; }

	const uint64_t kOnesStep4  = 0x1111111111111111ULL;
	const uint64_t kOnesStep8  = 0x0101010101010101ULL;
	const uint64_t kMSBsStep8  = 0x80ULL * kOnesStep8;

	uint64_t s = x;
	s = s - ((s & 0xA * kOnesStep4) >> 1);
	s = (s & 0x3 * kOnesStep4) + ((s >> 2) & 0x3 * kOnesStep4);
	s = (s + (s >> 4)) & 0xF * kOnesStep8;
	uint64_t byteSums = s * kOnesStep8;

	uint64_t kStep8 = k * kOnesStep8;
	uint64_t geqKStep8 = (((kStep8 | kMSBsStep8) - byteSums) & kMSBsStep8);
	uint64_t place = popcnt(geqKStep8) * 8;
	uint64_t byteRank = k - (((byteSums << 8) >> place) & (uint64_t)(0xFF));
	return place + kSelectInByte[((x >> place) & 0xFF) | (byteRank << 8)];
}

// Returns the position of the rank'th 1.  (rank = 0 returns the 1st 1)
// Returns 64 if there are fewer than rank+1 1s.
static inline uint64_t bitselect(uint64_t val, int rank) {
#ifdef __SSE4_2_
	uint64_t i = 1ULL << rank;
	asm("pdep %[val], %[mask], %[val]"
			: [val] "+r" (val)
			: [mask] "r" (i));
	asm("tzcnt %[bit], %[index]"
			: [index] "=r" (i)
			: [bit] "g" (val)
			: "cc");
	return i;
#endif
	return _select64(val, rank);
}

static inline uint64_t bitselectv(const uint64_t val, int ignore, int rank)
{
	return bitselect(val & ~BITMASK(ignore % 64), rank);
}

#if QF_BITS_PER_SLOT > 0
static inline qfblock * get_block(const QF *qf, uint64_t block_index)
{
	return &qf->blocks[block_index];
}
#else
static inline qfblock * get_block(const QF *qf, uint64_t block_index)
{
	return (qfblock *)(((char *)qf->blocks) + block_index * (sizeof(qfblock) +
						QF_SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8));
}
#endif

static inline int is_runend(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, runends, index) >> ((index % QF_SLOTS_PER_BLOCK) %
																								64)) & 1ULL;
}

static inline int is_occupied(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, occupieds, index) >> ((index % QF_SLOTS_PER_BLOCK) %
																									64)) & 1ULL;
}

#if QF_BITS_PER_SLOT == 8 || QF_BITS_PER_SLOT == 16 || QF_BITS_PER_SLOT == 32 || QF_BITS_PER_SLOT == 64

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	assert(index < qf->metadata->xnslots);
	return get_block(qf, index / QF_SLOTS_PER_BLOCK)->slots[index % QF_SLOTS_PER_BLOCK];
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->metadata->xnslots);
	get_block(qf, index / QF_SLOTS_PER_BLOCK)->slots[index % QF_SLOTS_PER_BLOCK] =
		value & BITMASK(qf->metadata->bits_per_slot);
}

#elif QF_BITS_PER_SLOT > 0

/* Little-endian code ....  Big-endian is TODO */

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	assert(index < qf->metadata->xnslots);
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 QF_SLOTS_PER_BLOCK)->slots[(index %
																																QF_SLOTS_PER_BLOCK)
																			 * QF_BITS_PER_SLOT / 8];
	return (uint64_t)(((*p) >> (((index % QF_SLOTS_PER_BLOCK) * QF_BITS_PER_SLOT) %
															8)) & BITMASK(QF_BITS_PER_SLOT));
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	assert(index < qf->metadata->xnslots);
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 QF_SLOTS_PER_BLOCK)->slots[(index %
																																QF_SLOTS_PER_BLOCK)
																			 * QF_BITS_PER_SLOT / 8];
	uint64_t t = *p;
	uint64_t mask = BITMASK(QF_BITS_PER_SLOT);
	uint64_t v = value;
	int shift = ((index % QF_SLOTS_PER_BLOCK) * QF_BITS_PER_SLOT) % 8;
	mask <<= shift;
	v <<= shift;
	t &= ~mask;
	t |= v;
	*p = t;
}

#else

/* Little-endian code ....  Big-endian is TODO */

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	assert(index < qf->metadata->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 QF_SLOTS_PER_BLOCK)->slots[(index %
																																QF_SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];
	// you cannot just do *p to get the value, undefined behavior
	uint64_t pvalue;
	memcpy(&pvalue,p,sizeof(pvalue));
	return (uint64_t)((pvalue >> (((index % QF_SLOTS_PER_BLOCK) *
															 qf->metadata->bits_per_slot) % 8)) &
										BITMASK(qf->metadata->bits_per_slot));
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->metadata->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 QF_SLOTS_PER_BLOCK)->slots[(index %
																																QF_SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];
	// This is undefined:
	//uint64_t t = *p;
	uint64_t t;
	memcpy(&t,p,sizeof(t));
	uint64_t mask = BITMASK(qf->metadata->bits_per_slot);
	uint64_t v = value;
	int shift = ((index % QF_SLOTS_PER_BLOCK) * qf->metadata->bits_per_slot) % 8;
	mask <<= shift;
	v <<= shift;
	t &= ~mask;
	t |= v;
	// this is undefined
	//*p = t;
	memcpy(p,&t,sizeof(t));
}

#endif

static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index);

static inline uint64_t block_offset(const QF *qf, uint64_t blockidx)
{
	/* If we have extended counters and a 16-bit (or larger) offset
		 field, then we can safely ignore the possibility of overflowing
		 that field. */
	if (sizeof(qf->blocks[0].offset) > 1 ||
			get_block(qf, blockidx)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
		return get_block(qf, blockidx)->offset;

	return run_end(qf, QF_SLOTS_PER_BLOCK * blockidx - 1) - QF_SLOTS_PER_BLOCK *
		blockidx + 1;
}

static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index)
{
	uint64_t bucket_block_index       = hash_bucket_index / QF_SLOTS_PER_BLOCK;
	uint64_t bucket_intrablock_offset = hash_bucket_index % QF_SLOTS_PER_BLOCK;
	uint64_t bucket_blocks_offset = block_offset(qf, bucket_block_index);

	uint64_t bucket_intrablock_rank   = bitrank(get_block(qf,
																				bucket_block_index)->occupieds[0],
																				bucket_intrablock_offset);

	if (bucket_intrablock_rank == 0) {
		if (bucket_blocks_offset <= bucket_intrablock_offset)
			return hash_bucket_index;
		else
			return QF_SLOTS_PER_BLOCK * bucket_block_index + bucket_blocks_offset - 1;
	}

	uint64_t runend_block_index  = bucket_block_index + bucket_blocks_offset /
		QF_SLOTS_PER_BLOCK;
	uint64_t runend_ignore_bits  = bucket_blocks_offset % QF_SLOTS_PER_BLOCK;
	uint64_t runend_rank         = bucket_intrablock_rank - 1;
	uint64_t runend_block_offset = bitselectv(get_block(qf,
																						runend_block_index)->runends[0],
																						runend_ignore_bits, runend_rank);
	if (runend_block_offset == QF_SLOTS_PER_BLOCK) {
		if (bucket_blocks_offset == 0 && bucket_intrablock_rank == 0) {
			/* The block begins in empty space, and this bucket is in that region of
			 * empty space */
			return hash_bucket_index;
		} else {
			do {
				runend_rank        -= popcntv(get_block(qf,
																								runend_block_index)->runends[0],
																			runend_ignore_bits);
				runend_block_index++;
				runend_ignore_bits  = 0;
				runend_block_offset = bitselectv(get_block(qf,
																									 runend_block_index)->runends[0],
																				 runend_ignore_bits, runend_rank);
			} while (runend_block_offset == QF_SLOTS_PER_BLOCK);
		}
	}

	uint64_t runend_index = QF_SLOTS_PER_BLOCK * runend_block_index +
		runend_block_offset;
	if (runend_index < hash_bucket_index)
		return hash_bucket_index;
	else
		return runend_index;
}

static inline int offset_lower_bound(const QF *qf, uint64_t slot_index)
{
	const qfblock * b = get_block(qf, slot_index / QF_SLOTS_PER_BLOCK);
	const uint64_t slot_offset = slot_index % QF_SLOTS_PER_BLOCK;
	const uint64_t boffset = b->offset;
	const uint64_t occupieds = b->occupieds[0] & BITMASK(slot_offset+1);
	assert(QF_SLOTS_PER_BLOCK == 64);
	if (boffset <= slot_offset) {
		const uint64_t runends = (b->runends[0] & BITMASK(slot_offset)) >> boffset;
		return popcnt(occupieds) - popcnt(runends);
	}
	return boffset - slot_offset + popcnt(occupieds);
}

static inline int qf_is_empty(const QF *qf, uint64_t slot_index)
{
	return offset_lower_bound(qf, slot_index) == 0;
}

static inline int might_be_empty(const QF *qf, uint64_t slot_index)
{
	return !is_occupied(qf, slot_index)
		&& !is_runend(qf, slot_index);
}

static inline int probably_is_empty(const QF *qf, uint64_t slot_index)
{
	return get_slot(qf, slot_index) == 0
		&& !is_occupied(qf, slot_index)
		&& !is_runend(qf, slot_index);
}

static inline uint64_t find_first_empty_slot(QF *qf, uint64_t from)
{
	do {
		int t = offset_lower_bound(qf, from);
		assert(t>=0);
		if (t == 0)
			break;
		from = from + t;
	} while(1);
	return from;
}

static inline uint64_t shift_into_b(const uint64_t a, const uint64_t b,
																		const int bstart, const int bend,
																		const int amount)
{
	const uint64_t a_component = bstart == 0 ? (a >> (64 - amount)) : 0;
	const uint64_t b_shifted_mask = BITMASK(bend - bstart) << bstart;
	const uint64_t b_shifted = ((b_shifted_mask & b) << amount) & b_shifted_mask;
	const uint64_t b_mask = ~b_shifted_mask;
	return a_component | b_shifted | (b & b_mask);
}

#if QF_BITS_PER_SLOT == 8 || QF_BITS_PER_SLOT == 16 || QF_BITS_PER_SLOT == 32 || QF_BITS_PER_SLOT == 64

static inline void shift_remainders(QF *qf, uint64_t start_index, uint64_t
																		empty_index)
{
	uint64_t start_block  = start_index / QF_SLOTS_PER_BLOCK;
	uint64_t start_offset = start_index % QF_SLOTS_PER_BLOCK;
	uint64_t empty_block  = empty_index / QF_SLOTS_PER_BLOCK;
	uint64_t empty_offset = empty_index % QF_SLOTS_PER_BLOCK;

	assert (start_index <= empty_index && empty_index < qf->metadata->xnslots);

	while (start_block < empty_block) {
		memmove(&get_block(qf, empty_block)->slots[1],
						&get_block(qf, empty_block)->slots[0],
						empty_offset * sizeof(qf->blocks[0].slots[0]));
		get_block(qf, empty_block)->slots[0] = get_block(qf,
																			empty_block-1)->slots[QF_SLOTS_PER_BLOCK-1];
		empty_block--;
		empty_offset = QF_SLOTS_PER_BLOCK-1;
	}

	memmove(&get_block(qf, empty_block)->slots[start_offset+1],
					&get_block(qf, empty_block)->slots[start_offset],
					(empty_offset - start_offset) * sizeof(qf->blocks[0].slots[0]));
}

#else

#define REMAINDER_WORD(qf, i) ((uint64_t *)&(get_block(qf, (i)/qf->metadata->bits_per_slot)->slots[8 * ((i) % qf->metadata->bits_per_slot)]))

static inline void shift_remainders(QF *qf, const uint64_t start_index, const
																		uint64_t empty_index)
{
	uint64_t last_word = (empty_index + 1) * qf->metadata->bits_per_slot / 64;
	const uint64_t first_word = start_index * qf->metadata->bits_per_slot / 64;
	int bend = ((empty_index + 1) * qf->metadata->bits_per_slot) % 64;
	const int bstart = (start_index * qf->metadata->bits_per_slot) % 64;

	while (last_word != first_word) {
		*REMAINDER_WORD(qf, last_word) = shift_into_b(*REMAINDER_WORD(qf, last_word-1),
																									*REMAINDER_WORD(qf, last_word),
																									0, bend, qf->metadata->bits_per_slot);
		last_word--;
		bend = 64;
	}
	*REMAINDER_WORD(qf, last_word) = shift_into_b(0, *REMAINDER_WORD(qf,
																																	 last_word),
																								bstart, bend,
																								qf->metadata->bits_per_slot);
}

#endif

/*
static inline void qf_dump_block(const QF *qf, uint64_t i)
{
	uint64_t j;

	printf("%-192d", get_block(qf, i)->offset);
	printf("\n");

	for (j = 0; j < QF_SLOTS_PER_BLOCK; j++)
		printf("%02lx ", j);
	printf("\n");

	for (j = 0; j < QF_SLOTS_PER_BLOCK; j++)
		printf(" %d ", (get_block(qf, i)->occupieds[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");

	for (j = 0; j < QF_SLOTS_PER_BLOCK; j++)
		printf(" %d ", (get_block(qf, i)->runends[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");

#if QF_BITS_PER_SLOT == 8 || QF_BITS_PER_SLOT == 16 || QF_BITS_PER_SLOT == 32
	for (j = 0; j < QF_SLOTS_PER_BLOCK; j++)
		printf("%02x ", get_block(qf, i)->slots[j]);
#elif QF_BITS_PER_SLOT == 64
	for (j = 0; j < QF_SLOTS_PER_BLOCK; j++)
		printf("%02lx ", get_block(qf, i)->slots[j]);
#else
	for (j = 0; j < QF_SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8; j++)
		printf("%02x ", get_block(qf, i)->slots[j]);
#endif

	printf("\n");

	printf("\n");
}

void qf_dump_metadata(const QF *qf) {
	printf("Slots: %lu Occupied: %lu Elements: %lu Distinct: %lu\n",
				 qf->metadata->nslots,
				 qf->metadata->noccupied_slots,
				 qf->metadata->nelts,
				 qf->metadata->ndistinct_elts);
	printf("Key_bits: %lu Value_bits: %lu Remainder_bits: %lu Bits_per_slot: %lu\n",
				 qf->metadata->key_bits,
				 qf->metadata->value_bits,
				 qf->metadata->key_remainder_bits,
				 qf->metadata->bits_per_slot);
}

void qf_dump(const QF *qf)
{
	uint64_t i;

	printf("%lu %lu %lu\n",
				 qf->metadata->nblocks,
				 qf->metadata->ndistinct_elts,
				 qf->metadata->nelts);

	for (i = 0; i < qf->metadata->nblocks; i++) {
		qf_dump_block(qf, i);
	}

}
*/

static inline void find_next_n_empty_slots(QF *qf, uint64_t from, uint64_t n,
																					 uint64_t *indices)
{
	while (n) {
		indices[--n] = find_first_empty_slot(qf, from);
		from = indices[n] + 1;
	}
}

static inline void shift_slots(QF *qf, int64_t first, uint64_t last, uint64_t
															 distance)
{
	int64_t i;
	if (distance == 1)
		shift_remainders(qf, first, last+1);
	else
		for (i = last; i >= first; i--)
			set_slot(qf, i + distance, get_slot(qf, i));
}

static inline void shift_runends(QF *qf, int64_t first, uint64_t last,
																 uint64_t distance)
{
	assert(last < qf->metadata->xnslots && distance < 64);
	uint64_t first_word = first / 64;
	uint64_t bstart = first % 64;
	uint64_t last_word = (last + distance + 1) / 64;
	uint64_t bend = (last + distance + 1) % 64;

	if (last_word != first_word) {
		METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(METADATA_WORD(qf, runends, 64*(last_word-1)),
																														METADATA_WORD(qf, runends, 64*last_word),
																														0, bend, distance);
		bend = 64;
		last_word--;
		while (last_word != first_word) {
			METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(METADATA_WORD(qf, runends, 64*(last_word-1)),
																															METADATA_WORD(qf, runends, 64*last_word),
																															0, bend, distance);
			last_word--;
		}
	}
	METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(0, METADATA_WORD(qf,
																																					 runends,
																																					 64*last_word),
																													bstart, bend, distance);

}

static inline bool insert_replace_slots_and_shift_remainders_and_runends_and_offsets(QF		*qf,
																																										 int		 operation,
																																										 uint64_t		 bucket_index,
																																										 uint64_t		 overwrite_index,
																																										 const uint64_t	*remainders,
																																										 uint64_t		 total_remainders,
																																										 uint64_t		 noverwrites)
{
	uint64_t empties[67];
	uint64_t i;
	int64_t j;
	int64_t ninserts = total_remainders - noverwrites;
	uint64_t insert_index = overwrite_index + noverwrites;

	if (ninserts > 0) {
		/* First, shift things to create n empty spaces where we need them. */
		find_next_n_empty_slots(qf, insert_index, ninserts, empties);
		if (empties[0] >= qf->metadata->xnslots) {
			return false;
		}
		for (j = 0; j < ninserts - 1; j++)
			shift_slots(qf, empties[j+1] + 1, empties[j] - 1, j + 1);
		shift_slots(qf, insert_index, empties[ninserts - 1] - 1, ninserts);

		for (j = 0; j < ninserts - 1; j++)
			shift_runends(qf, empties[j+1] + 1, empties[j] - 1, j + 1);
		shift_runends(qf, insert_index, empties[ninserts - 1] - 1, ninserts);

		for (i = noverwrites; i < total_remainders - 1; i++)
			METADATA_WORD(qf, runends, overwrite_index + i) &= ~(1ULL <<
																													 (((overwrite_index
																															+ i) %
																														 QF_SLOTS_PER_BLOCK)
																														% 64));

		switch (operation) {
			case 0: /* insert into empty bucket */
				assert (noverwrites == 0);
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |=
					1ULL << (((overwrite_index + total_remainders - 1) %
										QF_SLOTS_PER_BLOCK) % 64);
				break;
			case 1: /* append to bucket */
				METADATA_WORD(qf, runends, overwrite_index + noverwrites - 1)      &=
					~(1ULL << (((overwrite_index + noverwrites - 1) % QF_SLOTS_PER_BLOCK) %
										 64));
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |=
					1ULL << (((overwrite_index + total_remainders - 1) %
										QF_SLOTS_PER_BLOCK) % 64);
				break;
			case 2: /* insert into bucket */
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) &=
					~(1ULL << (((overwrite_index + total_remainders - 1) %
											QF_SLOTS_PER_BLOCK) % 64));
				break;
			default:
				fprintf(stderr, "Invalid operation %d\n", operation);
				abort();
		}

		uint64_t npreceding_empties = 0;
		for (i = bucket_index / QF_SLOTS_PER_BLOCK + 1; i <= empties[0]/QF_SLOTS_PER_BLOCK; i++) {
			while ((int64_t)npreceding_empties < ninserts &&
						 empties[ninserts - 1 - npreceding_empties]  / QF_SLOTS_PER_BLOCK < i)
				npreceding_empties++;

			if (get_block(qf, i)->offset + ninserts - npreceding_empties < BITMASK(8*sizeof(qf->blocks[0].offset)))
				get_block(qf, i)->offset += ninserts - npreceding_empties;
			else
				get_block(qf, i)->offset = (uint8_t) BITMASK(8*sizeof(qf->blocks[0].offset));
		}
	}

	for (i = 0; i < total_remainders; i++)
		set_slot(qf, overwrite_index + i, remainders[i]);

	modify_metadata(qf, &qf->metadata->noccupied_slots, ninserts);

	return true;
}

static inline int remove_replace_slots_and_shift_remainders_and_runends_and_offsets(QF		        *qf,
																																										 int		 operation,
																																										 uint64_t		 bucket_index,
																																										 uint64_t		 overwrite_index,
																																										 const uint64_t	*remainders,
																																										 uint64_t		 total_remainders,
																																										 uint64_t		 old_length)
{
	uint64_t i;

	// Update the slots
	for (i = 0; i < total_remainders; i++)
		set_slot(qf, overwrite_index + i, remainders[i]);

	// If this is the last thing in its run, then we may need to set a new runend bit
	if (is_runend(qf, overwrite_index + old_length - 1)) {
	  if (total_remainders > 0) {
	    // If we're not deleting this entry entirely, then it will still the last entry in this run
	    METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |= 1ULL << ((overwrite_index + total_remainders - 1) % 64);
	  } else if (overwrite_index > bucket_index &&
		     !is_runend(qf, overwrite_index - 1)) {
	    // If we're deleting this entry entirely, but it is not the first entry in this run,
	    // then set the preceding entry to be the runend
	    METADATA_WORD(qf, runends, overwrite_index - 1) |= 1ULL << ((overwrite_index - 1) % 64);
	  }
	}

	// shift slots back one run at a time
	uint64_t original_bucket = bucket_index;
	uint64_t current_bucket = bucket_index;
	uint64_t current_slot = overwrite_index + total_remainders;
	uint64_t current_distance = old_length - total_remainders;
	int ret_current_distance = current_distance;

	while (current_distance > 0) {
		if (is_runend(qf, current_slot + current_distance - 1)) {
			do {
				current_bucket++;
			} while (current_bucket < current_slot + current_distance &&
							 !is_occupied(qf, current_bucket));
		}

		if (current_bucket <= current_slot) {
			set_slot(qf, current_slot, get_slot(qf, current_slot + current_distance));
			if (is_runend(qf, current_slot) !=
					is_runend(qf, current_slot + current_distance))
				METADATA_WORD(qf, runends, current_slot) ^= 1ULL << (current_slot % 64);
			current_slot++;

		} else if (current_bucket <= current_slot + current_distance) {
			uint64_t i;
			for (i = current_slot; i < current_slot + current_distance; i++) {
				set_slot(qf, i, 0);
				METADATA_WORD(qf, runends, i) &= ~(1ULL << (i % 64));
			}

			current_distance = current_slot + current_distance - current_bucket;
			current_slot = current_bucket;
		} else {
			current_distance = 0;
		}
	}

	// reset the occupied bit of the hash bucket index if the hash is the
	// only item in the run and is removed completely.
	if (operation && !total_remainders)
		METADATA_WORD(qf, occupieds, bucket_index) &= ~(1ULL << (bucket_index % 64));

	// update the offset bits.
	// find the number of occupied slots in the original_bucket block.
	// Then find the runend slot corresponding to the last run in the
	// original_bucket block.
	// Update the offset of the block to which it belongs.
	uint64_t original_block = original_bucket / QF_SLOTS_PER_BLOCK;
	while (1 && old_length > total_remainders) {	// we only update offsets if we shift/delete anything
		int32_t last_occupieds_bit = bitscanreverse(get_block(qf, original_block)->occupieds[0]);
		// there is nothing in the block
		if (last_occupieds_bit == -1) {
			if (get_block(qf, original_block + 1)->offset == 0)
				break;
			get_block(qf, original_block + 1)->offset = 0;
		} else {
			uint64_t last_occupieds_hash_index = QF_SLOTS_PER_BLOCK * original_block + last_occupieds_bit;
			uint64_t runend_index = run_end(qf, last_occupieds_hash_index);
			// runend spans across the block
			// update the offset of the next block
			if (runend_index / QF_SLOTS_PER_BLOCK == original_block) { // if the run ends in the same block
				if (get_block(qf, original_block + 1)->offset == 0)
					break;
				get_block(qf, original_block + 1)->offset = 0;
			} else if (runend_index / QF_SLOTS_PER_BLOCK == original_block + 1) { // if the last run spans across one block
				if (get_block(qf, original_block + 1)->offset == (runend_index % QF_SLOTS_PER_BLOCK) + 1)
					break;
				get_block(qf, original_block + 1)->offset = (runend_index % QF_SLOTS_PER_BLOCK) + 1;
			} else { // if the last run spans across multiple blocks
				uint64_t i;
				for (i = original_block + 1; i < runend_index / QF_SLOTS_PER_BLOCK - 1; i++)
					get_block(qf, i)->offset = QF_SLOTS_PER_BLOCK;
				if (get_block(qf, runend_index / QF_SLOTS_PER_BLOCK)->offset == (runend_index % QF_SLOTS_PER_BLOCK) + 1)
					break;
				get_block(qf, runend_index / QF_SLOTS_PER_BLOCK)->offset = (runend_index % QF_SLOTS_PER_BLOCK) + 1;
			}
		}
		original_block++;
	}

	int num_slots_freed = old_length - total_remainders;
	modify_metadata(qf, &qf->metadata->noccupied_slots, -num_slots_freed);
	/*qf->metadata->noccupied_slots -= (old_length - total_remainders);*/
	if (!total_remainders) {
		modify_metadata(qf, &qf->metadata->ndistinct_elts, -1);
		/*qf->metadata->ndistinct_elts--;*/
	}

	return ret_current_distance;
}

/*****************************************************************************
 * Code that uses the above to implement a QF with keys and inline counters. *
 *****************************************************************************/

/*
	 Counter format:
	 0 xs:    <empty string>
	 1 x:     x
	 2 xs:    xx
	 3 0s:    000
	 >2 xs:   xbc...cx  for x != 0, b < x, c != 0, x
	 >3 0s:   0c...c00  for c != 0
	 */
static inline uint64_t *encode_counter(QF *qf, uint64_t remainder, uint64_t
																			 counter, uint64_t *slots)
{
	uint64_t digit = remainder;
	uint64_t base = (1ULL << qf->metadata->bits_per_slot) - 1;
	uint64_t *p = slots;

	if (counter == 0)
		return p;

	*--p = remainder;

	if (counter == 1)
		return p;

	if (counter == 2) {
		*--p = remainder;
		return p;
	}

	if (counter == 3 && remainder == 0) {
		*--p = remainder;
		*--p = remainder;
		return p;
	}

	if (counter == 3 && remainder > 0) {
		*--p = 0;
		*--p = remainder;
		return p;
	}

	if (remainder == 0)
		*--p = remainder;
	else
		base--;

	if (remainder)
		counter -= 3;
	else
		counter -= 4;
	do {
		digit = counter % base;
		digit++; /* Zero not allowed */
		if (remainder && digit >= remainder)
			digit++; /* Cannot overflow since digit is mod 2^r-2 */
		*--p = digit;
		counter /= base;
	} while (counter);

	if (remainder && digit >= remainder)
		*--p = 0;

	*--p = remainder;

	return p;
}

/* Returns the length of the encoding.
REQUIRES: index points to first slot of a counter. */
static inline uint64_t decode_counter(const QF *qf, uint64_t index, uint64_t
																			*remainder, uint64_t *count)
{
	uint64_t base;
	uint64_t rem;
	uint64_t cnt;
	uint64_t digit;
	uint64_t end;

	*remainder = rem = get_slot(qf, index);

	if (is_runend(qf, index)) { /* Entire run is "0" */
		*count = 1;
		return index;
	}

	digit = get_slot(qf, index + 1);

	if (is_runend(qf, index + 1)) {
		*count = digit == rem ? 2 : 1;
		return index + (digit == rem ? 1 : 0);
	}

	if (rem > 0 && digit >= rem) {
		*count = digit == rem ? 2 : 1;
		return index + (digit == rem ? 1 : 0);
	}

	if (rem > 0 && digit == 0 && get_slot(qf, index + 2) == rem) {
		*count = 3;
		return index + 2;
	}

	if (rem == 0 && digit == 0) {
		if (get_slot(qf, index + 2) == 0) {
			*count = 3;
			return index + 2;
		} else {
			*count = 2;
			return index + 1;
		}
	}

	cnt = 0;
	base = (1ULL << qf->metadata->bits_per_slot) - (rem ? 2 : 1);

	end = index + 1;
	while (digit != rem && !is_runend(qf, end)) {
		if (digit > rem)
			digit--;
		if (digit && rem)
			digit--;
		cnt = cnt * base + digit;

		end++;
		digit = get_slot(qf, end);
	}

	if (rem) {
		*count = cnt + 3;
		return end;
	}

	if (is_runend(qf, end) || get_slot(qf, end + 1) != 0) {
		*count = 1;
		return index;
	}

	*count = cnt + 4;
	return end + 1;
}

/* return the next slot which corresponds to a
 * different element
 * */
static inline uint64_t next_slot(QF *qf, uint64_t current)
{
	uint64_t rem = get_slot(qf, current);
	current++;

	while (get_slot(qf, current) == rem && current <= qf->metadata->nslots) {
		current++;
	}
	return current;
}

static inline int insert1(QF *qf, __uint128_t hash, uint8_t runtime_lock)
{
	int ret_distance = 0;
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot;
	uint64_t hash_bucket_block_offset = hash_bucket_index % QF_SLOTS_PER_BLOCK;

	if (GET_NO_LOCK(runtime_lock) != QF_NO_LOCK) {
		if (!qf_lock(qf, hash_bucket_index, /*small*/ true, runtime_lock))
			return QF_COULDNT_LOCK;
	}
	if (qf_is_empty(qf, hash_bucket_index) /* might_be_empty(qf, hash_bucket_index) && runend_index == hash_bucket_index */) {
		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		set_slot(qf, hash_bucket_index, hash_remainder);
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);

		ret_distance = 0;
		modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
		modify_metadata(qf, &qf->metadata->noccupied_slots, 1);
		modify_metadata(qf, &qf->metadata->nelts, 1);
	} else {
		uint64_t runend_index              = run_end(qf, hash_bucket_index);
		int operation = 0; /* Insert into empty bucket */
		uint64_t insert_index = runend_index + 1;
		uint64_t new_value = hash_remainder;

		/* printf("RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

		uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																	 hash_bucket_index
																																	 - 1) + 1;

		if (is_occupied(qf, hash_bucket_index)) {

			/* Find the counter for this remainder if it exists. */
			uint64_t current_remainder = get_slot(qf, runstart_index);
			uint64_t zero_terminator = runstart_index;

			/* The counter for 0 is special. */
			if (current_remainder == 0) {
				uint64_t t = runstart_index + 1;
				while (t < runend_index && get_slot(qf, t) != 0)
					t++;
				if (t < runend_index && get_slot(qf, t+1) == 0)
					zero_terminator = t+1; /* Three or more 0s */
				else if (runstart_index < runend_index && get_slot(qf, runstart_index
																													 + 1) == 0)
					zero_terminator = runstart_index + 1; /* Exactly two 0s */
				/* Otherwise, exactly one 0 (i.e. zero_terminator == runstart_index) */

				/* May read past end of run, but that's OK because loop below
					 can handle that */
				if (hash_remainder != 0) {
					runstart_index = zero_terminator + 1;
					current_remainder = get_slot(qf, runstart_index);
				}
			}

			/* Skip over counters for other remainders. */
			while (current_remainder < hash_remainder && runstart_index <=
						 runend_index) {
				/* If this remainder has an extended counter, skip over it. */
				if (runstart_index < runend_index &&
						get_slot(qf, runstart_index + 1) < current_remainder) {
					runstart_index = runstart_index + 2;
					while (runstart_index < runend_index &&
								 get_slot(qf, runstart_index) != current_remainder)
						runstart_index++;
					runstart_index++;

					/* This remainder has a simple counter. */
				} else {
					runstart_index++;
				}

				/* This may read past the end of the run, but the while loop
					 condition will prevent us from using the invalid result in
					 that case. */
				current_remainder = get_slot(qf, runstart_index);
			}

			/* If this is the first time we've inserted the new remainder,
				 and it is larger than any remainder in the run. */
			if (runstart_index > runend_index) {
				operation = 1;
				insert_index = runstart_index;
				new_value = hash_remainder;
				modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);

				/* This is the first time we're inserting this remainder, but
					 there are larger remainders already in the run. */
			} else if (current_remainder != hash_remainder) {
				operation = 2; /* Inserting */
				insert_index = runstart_index;
				new_value = hash_remainder;
				modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);

				/* Cases below here: we're incrementing the (simple or
					 extended) counter for this remainder. */

				/* If there's exactly one instance of this remainder. */
			} else if (runstart_index == runend_index ||
								 (hash_remainder > 0 && get_slot(qf, runstart_index + 1) >
									hash_remainder) ||
								 (hash_remainder == 0 && zero_terminator == runstart_index)) {
				operation = 2; /* Insert */
				insert_index = runstart_index;
				new_value = hash_remainder;

				/* If there are exactly two instances of this remainder. */
			} else if ((hash_remainder > 0 && get_slot(qf, runstart_index + 1) ==
									hash_remainder) ||
								 (hash_remainder == 0 && zero_terminator == runstart_index + 1)) {
				operation = 2; /* Insert */
				insert_index = runstart_index + 1;
				new_value = 0;

				/* Special case for three 0s */
			} else if (hash_remainder == 0 && zero_terminator == runstart_index + 2) {
				operation = 2; /* Insert */
				insert_index = runstart_index + 1;
				new_value = 1;

				/* There is an extended counter for this remainder. */
			} else {

				/* Move to the LSD of the counter. */
				insert_index = runstart_index + 1;
				while (get_slot(qf, insert_index+1) != hash_remainder)
					insert_index++;

				/* Increment the counter. */
				uint64_t digit, carry;
				do {
					carry = 0;
					digit = get_slot(qf, insert_index);
					// Convert a leading 0 (which is special) to a normal encoded digit
					if (digit == 0) {
						digit++;
						if (digit == current_remainder)
							digit++;
					}

					// Increment the digit
					digit = (digit + 1) & BITMASK(qf->metadata->bits_per_slot);

					// Ensure digit meets our encoding requirements
					if (digit == 0) {
						digit++;
						carry = 1;
					}
					if (digit == current_remainder)
						digit = (digit + 1) & BITMASK(qf->metadata->bits_per_slot);
					if (digit == 0) {
						digit++;
						carry = 1;
					}

					set_slot(qf, insert_index, digit);
					insert_index--;
				} while(insert_index > runstart_index && carry);

				/* If the counter needs to be expanded. */
				if (insert_index == runstart_index && (carry > 0 || (current_remainder
																														 != 0 && digit >=
																														 current_remainder)))
				{
					operation = 2; /* insert */
					insert_index = runstart_index + 1;
					if (!carry)						/* To prepend a 0 before the counter if the MSD is greater than the rem */
						new_value = 0;
					else if (carry) {			/* Increment the new value because we don't use 0 to encode counters */
						new_value = 2;
						/* If the rem is greater than or equal to the new_value then fail*/
						assert(new_value < current_remainder);
					}
				} else {
					operation = -1;
				}
			}
		} else {
			modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
		}

		if (operation >= 0) {
			uint64_t empty_slot_index = find_first_empty_slot(qf, runend_index+1);
			if (empty_slot_index >= qf->metadata->xnslots) {
				return QF_NO_SPACE;
			}
			shift_remainders(qf, insert_index, empty_slot_index);

			set_slot(qf, insert_index, new_value);
			ret_distance = insert_index - hash_bucket_index;

			shift_runends(qf, insert_index, empty_slot_index-1, 1);
			switch (operation) {
				case 0:
					METADATA_WORD(qf, runends, insert_index)   |= 1ULL << ((insert_index
																																	%
																																	QF_SLOTS_PER_BLOCK)
																																 % 64);
					break;
				case 1:
					METADATA_WORD(qf, runends, insert_index-1) &= ~(1ULL <<
																													(((insert_index-1) %
																														QF_SLOTS_PER_BLOCK) %
																													 64));
					METADATA_WORD(qf, runends, insert_index)   |= 1ULL << ((insert_index
																																	%
																																	QF_SLOTS_PER_BLOCK)
																																 % 64);
					break;
				case 2:
					METADATA_WORD(qf, runends, insert_index)   &= ~(1ULL <<
																													((insert_index %
																														QF_SLOTS_PER_BLOCK) %
																													 64));
					break;
				default:
					fprintf(stderr, "Invalid operation %d\n", operation);
					abort();
			}
			/*
			 * Increment the offset for each block between the hash bucket index
			 * and block of the empty slot
			 * */
			uint64_t i;
			for (i = hash_bucket_index / QF_SLOTS_PER_BLOCK + 1; i <=
					 empty_slot_index/QF_SLOTS_PER_BLOCK; i++) {
				if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
					get_block(qf, i)->offset++;
				assert(get_block(qf, i)->offset != 0);
			}
			modify_metadata(qf, &qf->metadata->noccupied_slots, 1);
		}
		modify_metadata(qf, &qf->metadata->nelts, 1);
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
	}

	if (GET_NO_LOCK(runtime_lock) != QF_NO_LOCK) {
		qf_unlock(qf, hash_bucket_index, /*small*/ true);
	}

	return ret_distance;
}

static inline int insert(QF *qf, __uint128_t hash, uint64_t count, uint8_t
												 runtime_lock)
{
	int ret_distance = 0;
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot;
	uint64_t hash_bucket_block_offset = hash_bucket_index % QF_SLOTS_PER_BLOCK;
	/*uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;*/

	if (GET_NO_LOCK(runtime_lock) != QF_NO_LOCK) {
		if (!qf_lock(qf, hash_bucket_index, /*small*/ false, runtime_lock))
			return QF_COULDNT_LOCK;
	}

	uint64_t runend_index             = run_end(qf, hash_bucket_index);

	/* Empty slot */
	if (might_be_empty(qf, hash_bucket_index) && runend_index ==
			hash_bucket_index) {
		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		set_slot(qf, hash_bucket_index, hash_remainder);
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);

		modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
		modify_metadata(qf, &qf->metadata->noccupied_slots, 1);
		modify_metadata(qf, &qf->metadata->nelts, 1);
		/* This trick will, I hope, keep the fast case fast. */
		if (count > 1) {
			insert(qf, hash, count - 1, QF_NO_LOCK);
		}
	} else { /* Non-empty slot */
		uint64_t new_values[67];
		int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																	hash_bucket_index
																																	- 1) + 1;

		bool ret;
		if (!is_occupied(qf, hash_bucket_index)) { /* Empty bucket, but its slot is occupied. */
			uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
			ret = insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																							0,
																																							hash_bucket_index,
																																							runstart_index,
																																							p,
																																							&new_values[67] - p,
																																							0);
			if (!ret)
				return QF_NO_SPACE;
			modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
			ret_distance = runstart_index - hash_bucket_index;
		} else { /* Non-empty bucket */

			uint64_t current_remainder, current_count, current_end;

			/* Find the counter for this remainder, if one exists. */
			current_end = decode_counter(qf, runstart_index, &current_remainder,
																	 &current_count);
			while (current_remainder < hash_remainder && !is_runend(qf, current_end)) {
				runstart_index = current_end + 1;
				current_end = decode_counter(qf, runstart_index, &current_remainder,
																		 &current_count);
			}

			/* If we reached the end of the run w/o finding a counter for this remainder,
				 then append a counter for this remainder to the run. */
			if (current_remainder < hash_remainder) {
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
				ret = insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																								1, /* Append to bucket */
																																								hash_bucket_index,
																																								current_end + 1,
																																								p,
																																								&new_values[67] - p,
																																								0);
				if (!ret)
					return QF_NO_SPACE;
				modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
				ret_distance = (current_end + 1) - hash_bucket_index;
				/* Found a counter for this remainder.  Add in the new count. */
			} else if (current_remainder == hash_remainder) {
				uint64_t *p = encode_counter(qf, hash_remainder, current_count + count, &new_values[67]);
				ret = insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																					is_runend(qf, current_end) ? 1 : 2,
																																					hash_bucket_index,
																																					runstart_index,
																																					p,
																																					&new_values[67] - p,
																																					current_end - runstart_index + 1);
			if (!ret)
				return QF_NO_SPACE;
			ret_distance = runstart_index - hash_bucket_index;
				/* No counter for this remainder, but there are larger
					 remainders, so we're not appending to the bucket. */
			} else {
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
				ret = insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																								2, /* Insert to bucket */
																																								hash_bucket_index,
																																								runstart_index,
																																								p,
																																								&new_values[67] - p,
																																								0);
				if (!ret)
					return QF_NO_SPACE;
				modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
			ret_distance = runstart_index - hash_bucket_index;
			}
		}
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL << (hash_bucket_block_offset % 64);

		modify_metadata(qf, &qf->metadata->nelts, count);
	}

	if (GET_NO_LOCK(runtime_lock) != QF_NO_LOCK) {
		qf_unlock(qf, hash_bucket_index, /*small*/ false);
	}

	return ret_distance;
}

inline static int _remove(QF *qf, __uint128_t hash, uint64_t count, uint8_t
													runtime_lock)
{
	int ret_numfreedslots = 0;
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot;
	uint64_t current_remainder, current_count, current_end;
	uint64_t new_values[67];

	if (GET_NO_LOCK(runtime_lock) != QF_NO_LOCK) {
		if (!qf_lock(qf, hash_bucket_index, /*small*/ false, runtime_lock))
			return -2;
	}

	/* Empty bucket */
	if (!is_occupied(qf, hash_bucket_index))
		return -1;

	uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf, hash_bucket_index - 1) + 1;
	uint64_t original_runstart_index = runstart_index;
	int only_item_in_the_run = 0;

	/*Find the counter for this remainder, if one exists.*/
	current_end = decode_counter(qf, runstart_index, &current_remainder, &current_count);
	while (current_remainder < hash_remainder && !is_runend(qf, current_end)) {
		runstart_index = current_end + 1;
		current_end = decode_counter(qf, runstart_index, &current_remainder, &current_count);
	}
	/* remainder not found in the given run */
	if (current_remainder != hash_remainder)
		return -1;

	if (original_runstart_index == runstart_index && is_runend(qf, current_end))
		only_item_in_the_run = 1;

	/* endode the new counter */
	uint64_t *p = encode_counter(qf, hash_remainder,
															 count > current_count ? 0 : current_count - count,
															 &new_values[67]);
	ret_numfreedslots = remove_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																		only_item_in_the_run,
																																		hash_bucket_index,
																																		runstart_index,
																																		p,
																																		&new_values[67] - p,
																																		current_end - runstart_index + 1);

	// update the nelements.
	modify_metadata(qf, &qf->metadata->nelts, -count);
	/*qf->metadata->nelts -= count;*/

	if (GET_NO_LOCK(runtime_lock) != QF_NO_LOCK) {
		qf_unlock(qf, hash_bucket_index, /*small*/ false);
	}

	return ret_numfreedslots;
}

/***********************************************************************
 * Code that uses the above to implement key-value-counter operations. *
 ***********************************************************************/

uint64_t qf_init(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t value_bits,
								 enum qf_hashmode hash, uint32_t seed, void* buffer, uint64_t
								 buffer_len)
{
	uint64_t num_slots, xnslots, nblocks;
	uint64_t key_remainder_bits, bits_per_slot;
	uint64_t size;
	uint64_t total_num_bytes;

	assert(popcnt(nslots) == 1); /* nslots must be a power of 2 */
	num_slots = nslots;
	xnslots = nslots + 10*sqrt((double)nslots);
	nblocks = (xnslots + QF_SLOTS_PER_BLOCK - 1) / QF_SLOTS_PER_BLOCK;
	key_remainder_bits = key_bits;
	while (nslots > 1) {
		assert(key_remainder_bits > 0);
		key_remainder_bits--;
		nslots >>= 1;
	}

	bits_per_slot = key_remainder_bits + value_bits;
	assert (QF_BITS_PER_SLOT == 0 || QF_BITS_PER_SLOT == qf->metadata->bits_per_slot);
	assert(bits_per_slot > 1);
#if QF_BITS_PER_SLOT == 8 || QF_BITS_PER_SLOT == 16 || QF_BITS_PER_SLOT == 32 || QF_BITS_PER_SLOT == 64
	size = nblocks * sizeof(qfblock);
#else
	size = nblocks * (sizeof(qfblock) + QF_SLOTS_PER_BLOCK * bits_per_slot / 8);
#endif

	total_num_bytes = sizeof(qfmetadata) + size;
	if (buffer == NULL || total_num_bytes > buffer_len)
		return total_num_bytes;
	memset(buffer, 0, total_num_bytes);
	qf->metadata = (qfmetadata *)(buffer);
	qf->blocks = (qfblock *)(qf->metadata + 1);

	qf->metadata->magic_endian_number = MAGIC_NUMBER;
	qf->metadata->auto_resize = 0;
	qf->metadata->hash_mode = hash;
	qf->metadata->total_size_in_bytes = size;
	qf->metadata->seed = seed;
	qf->metadata->nslots = num_slots;
	qf->metadata->xnslots = xnslots;
	qf->metadata->key_bits = key_bits;
	qf->metadata->value_bits = value_bits;
	qf->metadata->key_remainder_bits = key_remainder_bits;
	qf->metadata->bits_per_slot = bits_per_slot;

	qf->metadata->range = qf->metadata->nslots;
	qf->metadata->range <<= qf->metadata->key_remainder_bits;
	qf->metadata->nblocks = (qf->metadata->xnslots + QF_SLOTS_PER_BLOCK - 1) /
		QF_SLOTS_PER_BLOCK;
	qf->metadata->nelts = 0;
	qf->metadata->ndistinct_elts = 0;
	qf->metadata->noccupied_slots = 0;

	qf->runtimedata->num_locks = (qf->metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;
	qf->runtimedata->f_info.filepath = NULL;

	/* initialize all the locks to 0 */
	qf->runtimedata->metadata_lock = 0;
	qf->runtimedata->locks = (volatile int *)calloc(qf->runtimedata->num_locks,
																					sizeof(volatile int));
	if (qf->runtimedata->locks == NULL) {
		perror("Couldn't allocate memory for runtime locks.");
		exit(EXIT_FAILURE);
	}
#ifdef LOG_WAIT_TIME
	qf->runtimedata->wait_times = (wait_time_data*
																 )calloc(qf->runtimedata->num_locks+1,
																				 sizeof(wait_time_data));
	if (qf->runtimedata->wait_times == NULL) {
		perror("Couldn't allocate memory for runtime wait_times.");
		exit(EXIT_FAILURE);
	}
#endif

	return total_num_bytes;
}

uint64_t qf_use(QF* qf, void* buffer, uint64_t buffer_len)
{
	qf->metadata = (qfmetadata *)(buffer);
	if (qf->metadata->total_size_in_bytes + sizeof(qfmetadata) > buffer_len) {
		return qf->metadata->total_size_in_bytes + sizeof(qfmetadata);
	}
	qf->blocks = (qfblock *)(qf->metadata + 1);

	qf->runtimedata = (qfruntime *)calloc(sizeof(qfruntime), 1);
	if (qf->runtimedata == NULL) {
		perror("Couldn't allocate memory for runtime data.");
		exit(EXIT_FAILURE);
	}
	/* initialize all the locks to 0 */
	qf->runtimedata->metadata_lock = 0;
	qf->runtimedata->locks = (volatile int *)calloc(qf->runtimedata->num_locks,
																					sizeof(volatile int));
	if (qf->runtimedata->locks == NULL) {
		perror("Couldn't allocate memory for runtime locks.");
		exit(EXIT_FAILURE);
	}
#ifdef LOG_WAIT_TIME
	qf->runtimedata->wait_times = (wait_time_data*
																 )calloc(qf->runtimedata->num_locks+1,
																				 sizeof(wait_time_data));
	if (qf->runtimedata->wait_times == NULL) {
		perror("Couldn't allocate memory for runtime wait_times.");
		exit(EXIT_FAILURE);
	}
#endif

	return sizeof(qfmetadata) + qf->metadata->total_size_in_bytes;
}

void *qf_destroy(QF *qf)
{
	assert(qf->runtimedata->locks != NULL);
	free((void*)qf->runtimedata->locks);
	assert(qf->runtimedata != NULL);
	free(qf->runtimedata);

	return (void*)qf->metadata;
}

bool qf_malloc(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t
							 value_bits, enum qf_hashmode hash, uint32_t seed)
{
	uint64_t total_num_bytes = qf_init(qf, nslots, key_bits, value_bits,
																		 hash, seed, NULL, 0);

	void *buffer = malloc(total_num_bytes);
	if (buffer == NULL) {
		perror("Couldn't allocate memory for the CQF.");
		exit(EXIT_FAILURE);
	}

	qf->runtimedata = (qfruntime *)calloc(sizeof(qfruntime), 1);
	if (qf->runtimedata == NULL) {
		perror("Couldn't allocate memory for runtime data.");
		exit(EXIT_FAILURE);
	}

	uint64_t init_size = qf_init(qf, nslots, key_bits, value_bits, hash, seed,
															 buffer, total_num_bytes);

	if (init_size == total_num_bytes)
		return true;
	else
		return false;
}

bool qf_free(QF *qf)
{
	assert(qf->metadata != NULL);
	void *buffer = qf_destroy(qf);
	if (buffer != NULL) {
		free(buffer);
		return true;
	}

	return false;
}

void qf_copy(QF *dest, const QF *src)
{
	DEBUG_CQF("%s\n","Source CQF");
	DEBUG_DUMP(src);
	memcpy(dest->runtimedata, src->runtimedata, sizeof(qfruntime));
	memcpy(dest->metadata, src->metadata, sizeof(qfmetadata));
	memcpy(dest->blocks, src->blocks, src->metadata->total_size_in_bytes);
	DEBUG_CQF("%s\n","Destination CQF after copy.");
	DEBUG_DUMP(dest);
}

void qf_reset(QF *qf)
{
	qf->metadata->nelts = 0;
	qf->metadata->ndistinct_elts = 0;
	qf->metadata->noccupied_slots = 0;

#ifdef LOG_WAIT_TIME
	memset(qf->wait_times, 0,
				 (qf->runtimedata->num_locks+1)*sizeof(wait_time_data));
#endif
#if QF_BITS_PER_SLOT == 8 || QF_BITS_PER_SLOT == 16 || QF_BITS_PER_SLOT == 32 || QF_BITS_PER_SLOT == 64
	memset(qf->blocks, 0, qf->metadata->nblocks* sizeof(qfblock));
#else
	memset(qf->blocks, 0, qf->metadata->nblocks*(sizeof(qfblock) + QF_SLOTS_PER_BLOCK *
																		 qf->metadata->bits_per_slot / 8));
#endif
}

int64_t qf_resize_malloc(QF *qf, uint64_t nslots)
{
	QF new_qf;
	if (!qf_malloc(&new_qf, nslots, qf->metadata->key_bits,
								 qf->metadata->value_bits, qf->metadata->hash_mode,
								 qf->metadata->seed))
		return false;
	if (qf->metadata->auto_resize)
		qf_set_auto_resize(&new_qf, true);

	// copy keys from qf into new_qf
	QFi qfi;
	qf_iterator_from_position(qf, &qfi, 0);
	int64_t ret_numkeys = 0;
	do {
		uint64_t key, value, count;
		qfi_get_hash(&qfi, &key, &value, &count);
		qfi_next(&qfi);
		int ret = qf_insert(&new_qf, key, value, count, QF_NO_LOCK | QF_KEY_IS_HASH);
		if (ret < 0) {
			fprintf(stderr, "Failed to insert key: %" PRIx64 " into the new CQF.\n", key);
			return ret;
		}
		ret_numkeys++;
	} while(!qfi_end(&qfi));

	qf_free(qf);
	memcpy(qf, &new_qf, sizeof(QF));

	return ret_numkeys;
}

uint64_t qf_resize(QF* qf, uint64_t nslots, void* buffer, uint64_t buffer_len)
{
	QF new_qf;
	new_qf.runtimedata = (qfruntime *)calloc(sizeof(qfruntime), 1);
	if (new_qf.runtimedata == NULL) {
		perror("Couldn't allocate memory for runtime data.\n");
		exit(EXIT_FAILURE);
	}

	uint64_t init_size = qf_init(&new_qf, nslots, qf->metadata->key_bits,
															 qf->metadata->value_bits,
															 qf->metadata->hash_mode, qf->metadata->seed,
															 buffer, buffer_len);

	if (init_size > buffer_len)
		return init_size;

	if (qf->metadata->auto_resize)
		qf_set_auto_resize(&new_qf, true);

	// copy keys from qf into new_qf
	QFi qfi;
	qf_iterator_from_position(qf, &qfi, 0);
	do {
		uint64_t key, value, count;
		qfi_get_hash(&qfi, &key, &value, &count);
		qfi_next(&qfi);
		int ret = qf_insert(&new_qf, key, value, count, QF_NO_LOCK | QF_KEY_IS_HASH);
		if (ret < 0) {
			fprintf(stderr, "Failed to insert key: %" PRIx64 " into the new CQF.\n", key);
			abort();
		}
	} while(!qfi_end(&qfi));

	qf_free(qf);
	memcpy(qf, &new_qf, sizeof(QF));

	return init_size;
}

void qf_set_auto_resize(QF* qf, bool enabled)
{
	if (enabled)
		qf->metadata->auto_resize = 1;
	else
		qf->metadata->auto_resize = 0;
}

int qf_insert(QF *qf, uint64_t key, uint64_t value, uint64_t count, uint8_t
							flags)
{
	// We fill up the CQF up to 95% load factor.
	// This is a very conservative check.
	if (qf->metadata->noccupied_slots >= qf->metadata->nslots * 0.95) {
		if (qf->metadata->auto_resize) {
			fprintf(stdout, "Resizing the CQF.\n");
			qf_resize_malloc(qf, qf->metadata->nslots * 2);
		} else
			return QF_NO_SPACE;
	}
	if (count == 0)
		return 0;

	if (GET_KEY_HASH(flags) != QF_KEY_IS_HASH) {
		if (qf->metadata->hash_mode == QF_HASH_DEFAULT)
			key = MurmurHash64A(((void *)&key), sizeof(key),
													qf->metadata->seed) % qf->metadata->range;
		else if (qf->metadata->hash_mode == QF_HASH_INVERTIBLE)
			key = hash_64(key, BITMASK(qf->metadata->key_bits));
	}
	uint64_t hash = (key << qf->metadata->value_bits) | (value &
																											 BITMASK(qf->metadata->value_bits));
	int ret;
	if (count == 1)
		ret = insert1(qf, hash, flags);
	else
		ret = insert(qf, hash, count, flags);

	// check for fullness based on the distance from the home slot to the slot
	// in which the key is inserted
	if (ret > DISTANCE_FROM_HOME_SLOT_CUTOFF) {
		if (qf->metadata->auto_resize) {
			fprintf(stdout, "Resizing the CQF.\n");
			qf_resize_malloc(qf, qf->metadata->nslots * 2);
		} else {
			fprintf(stderr, "The CQF is filling up.\n");
		}
	}
	return ret;
}

int qf_set_count(QF *qf, uint64_t key, uint64_t value, uint64_t count, uint8_t
								 flags)
{
	if (count == 0)
		return 0;

	uint64_t cur_count = qf_count_key_value(qf, key, value, flags);
	int64_t delta = count - cur_count;

	int ret;
	if (delta == 0)
		ret = 0;
	else if (delta > 0)
		ret = qf_insert(qf, key, value, delta, flags);
	else
		ret = qf_remove(qf, key, value, labs(delta), flags);

	return ret;
}

int qf_remove(QF *qf, uint64_t key, uint64_t value, uint64_t count, uint8_t
							flags)
{
	if (count == 0)
		return true;

	if (GET_KEY_HASH(flags) != QF_KEY_IS_HASH) {
		if (qf->metadata->hash_mode == QF_HASH_DEFAULT)
			key = MurmurHash64A(((void *)&key), sizeof(key),
													qf->metadata->seed) % qf->metadata->range;
		else if (qf->metadata->hash_mode == QF_HASH_INVERTIBLE)
			key = hash_64(key, BITMASK(qf->metadata->key_bits));
	}
	uint64_t hash = (key << qf->metadata->value_bits) | (value &
																											 BITMASK(qf->metadata->value_bits));
	return _remove(qf, hash, count, flags);
}

int qf_delete_key_value(QF *qf, uint64_t key, uint64_t value, uint8_t flags)
{
	uint64_t count = qf_count_key_value(qf, key, value, flags);
	if (count == 0)
		return true;

	if (GET_KEY_HASH(flags) != QF_KEY_IS_HASH) {
		if (qf->metadata->hash_mode == QF_HASH_DEFAULT)
			key = MurmurHash64A(((void *)&key), sizeof(key),
													qf->metadata->seed) % qf->metadata->range;
		else if (qf->metadata->hash_mode == QF_HASH_INVERTIBLE)
			key = hash_64(key, BITMASK(qf->metadata->key_bits));
	}
	uint64_t hash = (key << qf->metadata->value_bits) | (value &
																											 BITMASK(qf->metadata->value_bits));
	return _remove(qf, hash, count, flags);
}

uint64_t qf_count_key_value(const QF *qf, uint64_t key, uint64_t value,
														uint8_t flags)
{
	if (GET_KEY_HASH(flags) != QF_KEY_IS_HASH) {
		if (qf->metadata->hash_mode == QF_HASH_DEFAULT)
			key = MurmurHash64A(((void *)&key), sizeof(key),
													qf->metadata->seed) % qf->metadata->range;
		else if (qf->metadata->hash_mode == QF_HASH_INVERTIBLE)
			key = hash_64(key, BITMASK(qf->metadata->key_bits));
	}
	uint64_t hash = (key << qf->metadata->value_bits) | (value &
																											 BITMASK(qf->metadata->value_bits));
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->bits_per_slot);
	int64_t hash_bucket_index = hash >> qf->metadata->bits_per_slot;

	if (!is_occupied(qf, hash_bucket_index))
		return 0;

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder)
			return current_count;
		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));

	return 0;
}

uint64_t qf_query(const QF *qf, uint64_t key, uint64_t *value, uint8_t flags)
{
	if (GET_KEY_HASH(flags) != QF_KEY_IS_HASH) {
		if (qf->metadata->hash_mode == QF_HASH_DEFAULT)
			key = MurmurHash64A(((void *)&key), sizeof(key),
													qf->metadata->seed) % qf->metadata->range;
		else if (qf->metadata->hash_mode == QF_HASH_INVERTIBLE)
			key = hash_64(key, BITMASK(qf->metadata->key_bits));
	}
	uint64_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->key_remainder_bits);
	int64_t hash_bucket_index = hash >> qf->metadata->key_remainder_bits;

	if (!is_occupied(qf, hash_bucket_index))
		return 0;

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder,
																 &current_count);
		*value = current_remainder & BITMASK(qf->metadata->value_bits);
		current_remainder = current_remainder >> qf->metadata->value_bits;
		if (current_remainder == hash_remainder) {
			return current_count;
		}
		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));

	return 0;
}

int64_t qf_get_unique_index(const QF *qf, uint64_t key, uint64_t value,
														uint8_t flags)
{
	if (GET_KEY_HASH(flags) != QF_KEY_IS_HASH) {
		if (qf->metadata->hash_mode == QF_HASH_DEFAULT)
			key = MurmurHash64A(((void *)&key), sizeof(key),
													qf->metadata->seed) % qf->metadata->range;
		else if (qf->metadata->hash_mode == QF_HASH_INVERTIBLE)
			key = hash_64(key, BITMASK(qf->metadata->key_bits));
	}
	uint64_t hash = (key << qf->metadata->value_bits) | (value &
																											 BITMASK(qf->metadata->value_bits));
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->bits_per_slot);
	int64_t hash_bucket_index = hash >> qf->metadata->bits_per_slot;

	if (!is_occupied(qf, hash_bucket_index))
		return QF_DOESNT_EXIST;

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder)
			return runstart_index;

		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));

	return QF_DOESNT_EXIST;
}

enum qf_hashmode qf_get_hashmode(const QF *qf) {
	return qf->metadata->hash_mode;
}
uint64_t qf_get_hash_seed(const QF *qf) {
	return qf->metadata->seed;
}
__uint128_t qf_get_hash_range(const QF *qf) {
	return qf->metadata->range;
}

bool qf_is_auto_resize_enabled(const QF *qf) {
	if (qf->metadata->auto_resize == 1)
		return true;
	return false;
}
uint64_t qf_get_total_size_in_bytes(const QF *qf) {
	return qf->metadata->total_size_in_bytes;
}
uint64_t qf_get_nslots(const QF *qf) {
	return qf->metadata->nslots;
}
uint64_t qf_get_num_occupied_slots(const QF *qf) {
	return qf->metadata->noccupied_slots;
}

uint64_t qf_get_num_key_bits(const QF *qf) {
	return qf->metadata->key_bits;
}
uint64_t qf_get_num_value_bits(const QF *qf) {
	return qf->metadata->value_bits;
}
uint64_t qf_get_num_key_remainder_bits(const QF *qf) {
	return qf->metadata->key_remainder_bits;
}
uint64_t qf_get_bits_per_slot(const QF *qf) {
	return qf->metadata->bits_per_slot;
}

uint64_t qf_get_sum_of_counts(const QF *qf) {
	return qf->metadata->nelts;
}
uint64_t qf_get_num_distinct_key_value_pairs(const QF *qf) {
	return qf->metadata->ndistinct_elts;
}

/* initialize the iterator at the run corresponding
 * to the position index
 */
int64_t qf_iterator_from_position(const QF *qf, QFi *qfi, uint64_t position)
{
	if (position == 0xffffffffffffffff) {
		qfi->current = 0xffffffffffffffff;
		qfi->qf = qf;
		return QFI_INVALID;
	}
	assert(position < qf->metadata->nslots);
	if (!is_occupied(qf, position)) {
		uint64_t block_index = position;
		uint64_t idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
		if (idx == 64) {
			while(idx == 64 && block_index < qf->metadata->nblocks) {
				block_index++;
				idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
			}
		}
		position = block_index * QF_SLOTS_PER_BLOCK + idx;
	}

	qfi->qf = qf;
	qfi->num_clusters = 0;
	qfi->run = position;
	qfi->current = position == 0 ? 0 : run_end(qfi->qf, position-1) + 1;
	if (qfi->current < position)
		qfi->current = position;

#ifdef LOG_CLUSTER_LENGTH
	qfi->c_info = (cluster_data* )calloc(qf->metadata->nslots/32,
																			 sizeof(cluster_data));
	if (qfi->c_info == NULL) {
		perror("Couldn't allocate memory for c_info.");
		exit(EXIT_FAILURE);
	}
	qfi->cur_start_index = position;
	qfi->cur_length = 1;
#endif

	if (qfi->current >= qf->metadata->nslots)
		return QFI_INVALID;
	return qfi->current;
}

int64_t qf_iterator_key_value(const QF *qf, QFi *qfi, uint64_t key, uint64_t
															value, uint8_t flags)
{
	if (key >= qf->metadata->range) {
		qfi->current = 0xffffffffffffffff;
		qfi->qf = qf;
		return QFI_INVALID;
	}

	qfi->qf = qf;
	qfi->num_clusters = 0;

	if (GET_KEY_HASH(flags) != QF_KEY_IS_HASH) {
		if (qf->metadata->hash_mode == QF_HASH_DEFAULT)
			key = MurmurHash64A(((void *)&key), sizeof(key),
													qf->metadata->seed) % qf->metadata->range;
		else if (qf->metadata->hash_mode == QF_HASH_INVERTIBLE)
			key = hash_64(key, BITMASK(qf->metadata->key_bits));
	}
	uint64_t hash = (key << qf->metadata->value_bits) | (value &
																											 BITMASK(qf->metadata->value_bits));

	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index = hash >> qf->metadata->bits_per_slot;
	bool flag = false;

	// If a run starts at "position" move the iterator to point it to the
	// smallest key greater than or equal to "hash".
	if (is_occupied(qf, hash_bucket_index)) {
		uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																	 hash_bucket_index-1)
			+ 1;
		if (runstart_index < hash_bucket_index)
			runstart_index = hash_bucket_index;
		uint64_t current_remainder, current_count, current_end;
		do {
			current_end = decode_counter(qf, runstart_index, &current_remainder,
																	 &current_count);
			if (current_remainder >= hash_remainder) {
				flag = true;
				break;
			}
			runstart_index = current_end + 1;
		} while (!is_runend(qf, current_end));
		// found "hash" or smallest key greater than "hash" in this run.
		if (flag) {
			qfi->run = hash_bucket_index;
			qfi->current = runstart_index;
		}
	}
	// If a run doesn't start at "position" or the largest key in the run
	// starting at "position" is smaller than "hash" then find the start of the
	// next run.
	if (!is_occupied(qf, hash_bucket_index) || !flag) {
		uint64_t position = hash_bucket_index;
		assert(position < qf->metadata->nslots);
		uint64_t block_index = position / QF_SLOTS_PER_BLOCK;
		uint64_t idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
		if (idx == 64) {
			while(idx == 64 && block_index < qf->metadata->nblocks) {
				block_index++;
				idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
			}
		}
		position = block_index * QF_SLOTS_PER_BLOCK + idx;
		qfi->run = position;
		qfi->current = position == 0 ? 0 : run_end(qfi->qf, position-1) + 1;
		if (qfi->current < position)
			qfi->current = position;
	}

	if (qfi->current >= qf->metadata->nslots)
		return QFI_INVALID;
	return qfi->current;
}

static int qfi_get(const QFi *qfi, uint64_t *key, uint64_t *value, uint64_t
									 *count)
{
	if (qfi_end(qfi))
		return QFI_INVALID;

	uint64_t current_remainder, current_count;
	decode_counter(qfi->qf, qfi->current, &current_remainder, &current_count);

	*value = current_remainder & BITMASK(qfi->qf->metadata->value_bits);
	current_remainder = current_remainder >> qfi->qf->metadata->value_bits;
	*key = (qfi->run << qfi->qf->metadata->key_remainder_bits) | current_remainder;
	*count = current_count;

	return 0;
}

int qfi_get_key(const QFi *qfi, uint64_t *key, uint64_t *value, uint64_t
								*count)
{
	*key = *value = *count = 0;
	int ret = qfi_get(qfi, key, value, count);
	if (ret == 0) {
		if (qfi->qf->metadata->hash_mode == QF_HASH_DEFAULT) {
			*key = 0; *value = 0; *count = 0;
			return QF_INVALID;
		} else if (qfi->qf->metadata->hash_mode == QF_HASH_INVERTIBLE)
			*key = hash_64i(*key, BITMASK(qfi->qf->metadata->key_bits));
	}

	return ret;
}

int qfi_get_hash(const QFi *qfi, uint64_t *key, uint64_t *value, uint64_t
								 *count)
{
	*key = *value = *count = 0;
	return qfi_get(qfi, key, value, count);
}

int qfi_next(QFi *qfi)
{
	if (qfi_end(qfi))
		return QFI_INVALID;
	else {
		/* move to the end of the current counter*/
		uint64_t current_remainder, current_count;
		qfi->current = decode_counter(qfi->qf, qfi->current, &current_remainder,
																	&current_count);

		if (!is_runend(qfi->qf, qfi->current)) {
			qfi->current++;
#ifdef LOG_CLUSTER_LENGTH
			qfi->cur_length++;
#endif
			if (qfi_end(qfi))
				return QFI_INVALID;
			return 0;
		} else {
#ifdef LOG_CLUSTER_LENGTH
			/* save to check if the new current is the new cluster. */
			uint64_t old_current = qfi->current;
#endif
			uint64_t block_index = qfi->run / QF_SLOTS_PER_BLOCK;
			uint64_t rank = bitrank(get_block(qfi->qf, block_index)->occupieds[0],
															qfi->run % QF_SLOTS_PER_BLOCK);
			uint64_t next_run = bitselect(get_block(qfi->qf,
																							block_index)->occupieds[0],
																		rank);
			if (next_run == 64) {
				rank = 0;
				while (next_run == 64 && block_index < qfi->qf->metadata->nblocks) {
					block_index++;
					next_run = bitselect(get_block(qfi->qf, block_index)->occupieds[0],
															 rank);
				}
			}
			if (block_index == qfi->qf->metadata->nblocks) {
				/* set the index values to max. */
				qfi->run = qfi->current = qfi->qf->metadata->xnslots;
				return QFI_INVALID;
			}
			qfi->run = block_index * QF_SLOTS_PER_BLOCK + next_run;
			qfi->current++;
			if (qfi->current < qfi->run)
				qfi->current = qfi->run;
#ifdef LOG_CLUSTER_LENGTH
			if (qfi->current > old_current + 1) { /* new cluster. */
				if (qfi->cur_length > 10) {
					qfi->c_info[qfi->num_clusters].start_index = qfi->cur_start_index;
					qfi->c_info[qfi->num_clusters].length = qfi->cur_length;
					qfi->num_clusters++;
				}
				qfi->cur_start_index = qfi->run;
				qfi->cur_length = 1;
			} else {
				qfi->cur_length++;
			}
#endif
			return 0;
		}
	}
}

bool qfi_end(const QFi *qfi)
{
	if (qfi->current >= qfi->qf->metadata->xnslots /*&& is_runend(qfi->qf, qfi->current)*/)
		return true;
	return false;
}

/*
 * Merge qfa and qfb into qfc
 */
/*
 * iterate over both qf (qfa and qfb)
 * simultaneously
 * for each index i
 * min(get_value(qfa, ia) < get_value(qfb, ib))
 * insert(min, ic)
 * increment either ia or ib, whichever is minimum.
 */
void qf_merge(const QF *qfa, const QF *qfb, QF *qfc)
{
	QFi qfia, qfib;
	qf_iterator_from_position(qfa, &qfia, 0);
	qf_iterator_from_position(qfb, &qfib, 0);

	if (qfa->metadata->hash_mode != qfc->metadata->hash_mode &&
			qfa->metadata->seed != qfc->metadata->seed &&
			qfb->metadata->hash_mode  != qfc->metadata->hash_mode &&
			qfb->metadata->seed  != qfc->metadata->seed) {
		fprintf(stderr, "Output QF and input QFs do not have the same hash mode or seed.\n");
		exit(1);
	}

	uint64_t keya, valuea, counta, keyb, valueb, countb;
	qfi_get_hash(&qfia, &keya, &valuea, &counta);
	qfi_get_hash(&qfib, &keyb, &valueb, &countb);
	do {
		if (keya < keyb) {
			qf_insert(qfc, keya, valuea, counta, QF_NO_LOCK | QF_KEY_IS_HASH);
			qfi_next(&qfia);
			qfi_get_hash(&qfia, &keya, &valuea, &counta);
		}
		else {
			qf_insert(qfc, keyb, valueb, countb, QF_NO_LOCK | QF_KEY_IS_HASH);
			qfi_next(&qfib);
			qfi_get_hash(&qfib, &keyb, &valueb, &countb);
		}
	} while(!qfi_end(&qfia) && !qfi_end(&qfib));

	if (!qfi_end(&qfia)) {
		do {
			qfi_get_hash(&qfia, &keya, &valuea, &counta);
			qf_insert(qfc, keya, valuea, counta, QF_NO_LOCK | QF_KEY_IS_HASH);
		} while(!qfi_next(&qfia));
	}
	if (!qfi_end(&qfib)) {
		do {
			qfi_get_hash(&qfib, &keyb, &valueb, &countb);
			qf_insert(qfc, keyb, valueb, countb, QF_NO_LOCK | QF_KEY_IS_HASH);
		} while(!qfi_next(&qfib));
	}
}

/*
 * Merge an array of qfs into the resultant QF
 */
void qf_multi_merge(const QF *qf_arr[], int nqf, QF *qfr)
{
	int i;
	QFi qfi_arr[nqf];
	int smallest_idx = 0;
	uint64_t smallest_key = UINT64_MAX;
	for (i=0; i<nqf; i++) {
		if (qf_arr[i]->metadata->hash_mode != qfr->metadata->hash_mode &&
				qf_arr[i]->metadata->seed != qfr->metadata->seed) {
			fprintf(stderr, "Output QF and input QFs do not have the same hash mode or seed.\n");
			exit(1);
		}
		qf_iterator_from_position(qf_arr[i], &qfi_arr[i], 0);
	}

	DEBUG_CQF("Merging %d CQFs\n", nqf);
	for (i=0; i<nqf; i++) {
		DEBUG_CQF("CQF %d\n", i);
		DEBUG_DUMP(qf_arr[i]);
	}

	while (nqf > 1) {
		uint64_t keys[nqf];
		uint64_t values[nqf];
		uint64_t counts[nqf];
		for (i=0; i<nqf; i++)
			qfi_get_hash(&qfi_arr[i], &keys[i], &values[i], &counts[i]);

		do {
			smallest_key = UINT64_MAX;
			for (i=0; i<nqf; i++) {
				if (keys[i] < smallest_key) {
					smallest_key = keys[i]; smallest_idx = i;
				}
			}
			qf_insert(qfr, keys[smallest_idx], values[smallest_idx],
								counts[smallest_idx], QF_NO_LOCK | QF_KEY_IS_HASH);
			qfi_next(&qfi_arr[smallest_idx]);
			qfi_get_hash(&qfi_arr[smallest_idx], &keys[smallest_idx],
									 &values[smallest_idx],
							&counts[smallest_idx]);
		} while(!qfi_end(&qfi_arr[smallest_idx]));

		/* remove the qf that is exhausted from the array */
		if (smallest_idx < nqf-1)
			memmove(&qfi_arr[smallest_idx], &qfi_arr[smallest_idx+1],
							(nqf-smallest_idx-1)*sizeof(qfi_arr[0]));
		nqf--;
	}
	if (!qfi_end(&qfi_arr[0])) {
		uint64_t iters = 0;
		do {
			uint64_t key, value, count;
			qfi_get_hash(&qfi_arr[0], &key, &value, &count);
			qf_insert(qfr, key, value, count, QF_NO_LOCK | QF_KEY_IS_HASH);
			qfi_next(&qfi_arr[0]);
			iters++;
		} while(!qfi_end(&qfi_arr[0]));
		DEBUG_CQF("Num of iterations: %" PRIx64 "\n", iters);
	}

	DEBUG_CQF("%s", "Final CQF after merging.\n");
	DEBUG_DUMP(qfr);

	return;
}

/* find cosine similarity between two QFs. */
uint64_t qf_inner_product(const QF *qfa, const QF *qfb)
{
	uint64_t acc = 0;
	QFi qfi;
	const QF *qf_mem, *qf_disk;

	if (qfa->metadata->hash_mode != qfb->metadata->hash_mode &&
			qfa->metadata->seed != qfb->metadata->seed) {
		fprintf(stderr, "Input QFs do not have the same hash mode or seed.\n");
		exit(1);
	}

	// create the iterator on the larger QF.
	if (qfa->metadata->total_size_in_bytes > qfb->metadata->total_size_in_bytes)
	{
		qf_mem = qfb;
		qf_disk = qfa;
	} else {
		qf_mem = qfa;
		qf_disk = qfb;
	}

	qf_iterator_from_position(qf_disk, &qfi, 0);
	do {
		uint64_t key = 0, value = 0, count = 0;
		uint64_t count_mem;
		qfi_get_hash(&qfi, &key, &value, &count);
		if ((count_mem = qf_count_key_value(qf_mem, key, 0, QF_KEY_IS_HASH)) > 0) {
			acc += count*count_mem;
		}
	} while (!qfi_next(&qfi));

	return acc;
}

/* find cosine similarity between two QFs. */
void qf_intersect(const QF *qfa, const QF *qfb, QF *qfr)
{
	QFi qfi;
	const QF *qf_mem, *qf_disk;

	if (qfa->metadata->hash_mode != qfr->metadata->hash_mode &&
			qfa->metadata->seed != qfr->metadata->seed &&
			qfb->metadata->hash_mode  != qfr->metadata->hash_mode &&
			qfb->metadata->seed  != qfr->metadata->seed) {
		fprintf(stderr, "Output QF and input QFs do not have the same hash mode or seed.\n");
		exit(1);
	}

	// create the iterator on the larger QF.
	if (qfa->metadata->total_size_in_bytes > qfb->metadata->total_size_in_bytes)
	{
		qf_mem = qfb;
		qf_disk = qfa;
	} else {
		qf_mem = qfa;
		qf_disk = qfb;
	}

	qf_iterator_from_position(qf_disk, &qfi, 0);
	do {
		uint64_t key = 0, value = 0, count = 0;
		qfi_get_hash(&qfi, &key, &value, &count);
		if (qf_count_key_value(qf_mem, key, 0, QF_KEY_IS_HASH) > 0)
			qf_insert(qfr, key, value, count, QF_NO_LOCK | QF_KEY_IS_HASH);
	} while (!qfi_next(&qfi));
}

/* magnitude of a QF. */
uint64_t qf_magnitude(const QF *qf)
{
	return sqrt(qf_inner_product(qf, qf));
}

