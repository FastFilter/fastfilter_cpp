/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *
 * ============================================================================
 */

#ifndef _GQF_H_
#define _GQF_H_

#include <inttypes.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

	typedef struct quotient_filter quotient_filter;
	typedef quotient_filter QF;

	/* CQFs support three hashing modes:

		 - DEFAULT uses a hash that may introduce false positives, but
       this can be useful when inserting large keys that need to be
       hashed down to a small fingerprint.  With this type of hash,
       you can iterate over the hash values of all the keys in the
       CQF, but you cannot iterate over the keys themselves.

		 - INVERTIBLE has no false positives, but the size of the hash
       output must be the same as the size of the hash input,
       e.g. 17-bit keys hashed to 17-bit outputs.  So this mode is
       generally only useful when storing small keys in the CQF.  With
       this hashing mode, you can use iterators to enumerate both all
       the hashes in the CQF, or all the keys.

		 - NONE, for when you've done the hashing yourself.  WARNING: the
		   CQF can exhibit very bad performance if you insert a skewed
			 distribution of intputs.
	*/
	
	enum qf_hashmode {
		QF_HASH_DEFAULT,
		QF_HASH_INVERTIBLE,
		QF_HASH_NONE
	};

	/* The CQF supports concurrent insertions and queries.  Only the
		 portion of the CQF being examined or modified is locked, so it
		 supports high throughput even with many threads.

		 The CQF operations support 3 locking modes:

		 - NO_LOCK: for single-threaded applications or applications
       that do their own concurrency management.

		 - WAIT_FOR_LOCK: Spin until you get the lock, then do the query
       or update.

		 - TRY_ONCE_LOCK: If you can't grab the lock on the first try,
       return with an error code.
	*/
#define QF_NO_LOCK (0x01)
#define QF_TRY_ONCE_LOCK (0x02)
#define QF_WAIT_FOR_LOCK (0x04)

	/* It is sometimes useful to insert a key that has already been
		 hashed. */
#define QF_KEY_IS_HASH (0x08)

	/******************************************
		 The CQF defines low-level constructor and destructor operations
		 that are designed to enable the application to manage the memory
		 used by the CQF. 
	*******************************************/
	
	/*
	 * Create an empty CQF in "buffer".  If there is not enough space at
	 * buffer then it will return the total size needed in bytes to
	 * initialize the CQF.  This function takes ownership of buffer.
	 */
	uint64_t qf_init(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t
									 value_bits, enum qf_hashmode hash, uint32_t seed, void*
									 buffer, uint64_t buffer_len);

	/* Create a CQF in "buffer". Note that this does not initialize the
	 contents of bufferss Use this function if you have read a CQF, e.g.
	 off of disk or network, and want to begin using that stream of
	 bytes as a CQF. The CQF takes ownership of buffer.  */
	uint64_t qf_use(QF* qf, void* buffer, uint64_t buffer_len);

	/* Destroy this CQF.  Returns a pointer to the memory that the CQF was
		 using (i.e. passed into qf_init or qf_use) so that the application
		 can release that memory. */
	void *qf_destroy(QF *qf);

	/* Allocate a new CQF using "nslots" at "buffer" and copy elements from "qf"
	 * into it. 
	 * If there is not enough space at buffer then it will return the total size
	 * needed in bytes to initialize the new CQF.
	 * */
	uint64_t qf_resize(QF* qf, uint64_t nslots, void* buffer, uint64_t
										 buffer_len);

	/***********************************
    The following convenience functions create and destroy CQFs by
		using malloc/free to obtain and release the memory for the CQF. 
	************************************/
	
	/* Initialize the CQF and allocate memory for the CQF. */
	bool qf_malloc(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t
								 value_bits, enum qf_hashmode hash, uint32_t seed);

	bool qf_free(QF *qf);

	/* Resize the QF to the specified number of slots.  Uses malloc() to
	 * obtain the new memory, and calls free() on the old memory.
	 * Return value:
	 *    >= 0: number of keys copied during resizing.
	 * */
	int64_t qf_resize_malloc(QF *qf, uint64_t nslots);

	/* Turn on automatic resizing.  Resizing is performed by calling
		 qf_resize_malloc, so the CQF must meet the requirements of that
		 function. */
	void qf_set_auto_resize(QF* qf, bool enabled);

	/***********************************
   Functions for modifying the CQF.
	***********************************/

#define QF_NO_SPACE (-1)
#define QF_COULDNT_LOCK (-2)
#define QF_DOESNT_EXIST (-3)
	
	/* Increment the counter for this key/value pair by count. 
	 * Return value:
	 *    >= 0: distance from the home slot to the slot in which the key is
	 *          inserted (or 0 if count == 0).
	 *    == QF_NO_SPACE: the CQF has reached capacity.
	 *    == QF_COULDNT_LOCK: TRY_ONCE_LOCK has failed to acquire the lock.
	 */
	int qf_insert(QF *qf, uint64_t key, uint64_t value, uint64_t count, uint8_t
								flags);

	/* Set the counter for this key/value pair to count. 
	 Return value: Same as qf_insert. 
	 Returns 0 if new count is equal to old count.
	*/
	int qf_set_count(QF *qf, uint64_t key, uint64_t value, uint64_t count,
									 uint8_t flags);

	/* Remove up to count instances of this key/value combination.
	 * If the CQF contains <= count instances, then they will all be 
	 * removed, which is not an error.
	 * Return value:
	 *    >=  0: number of slots freed.
	 *    == QF_DOESNT_EXIST: Specified item did not exist.
	 *    == QF_COULDNT_LOCK: TRY_ONCE_LOCK has failed to acquire the lock.
	 */
	int qf_remove(QF *qf, uint64_t key, uint64_t value, uint64_t count, uint8_t
								flags);

	/* Remove all instances of this key/value pair. */
	int qf_delete_key_value(QF *qf, uint64_t key, uint64_t value, uint8_t flags);

	/* Remove all instances of this key. */
	/* NOT IMPLEMENTED YET. */
	//void qf_delete_key(QF *qf, uint64_t key);

	/* Replace the association (key, oldvalue, count) with the association
		 (key, newvalue, count). If there is already an association (key,
		 newvalue, count'), then the two associations will be merged and
		 their counters will be summed, resulting in association (key,
		 newvalue, count' + count). */
	/* NOT IMPLEMENTED YET. */
	//void qf_replace(QF *qf, uint64_t key, uint64_t oldvalue, uint64_t newvalue);

	/****************************************
   Query functions
	****************************************/
	
	/* Lookup the value associated with key.  Returns the count of that
		 key/value pair in the QF.  If it returns 0, then, the key is not
		 present in the QF. Only returns the first value associated with key
		 in the QF.  If you want to see others, use an iterator. 
		 May return QF_COULDNT_LOCK if called with QF_TRY_LOCK.  */
	uint64_t qf_query(const QF *qf, uint64_t key, uint64_t *value, uint8_t
										flags);

	/* Return the number of times key has been inserted, with any value,
		 into qf. */
	/* NOT IMPLEMENTED YET. */
	//uint64_t qf_count_key(const QF *qf, uint64_t key);

	/* Return the number of times key has been inserted, with the given
		 value, into qf.
		 May return QF_COULDNT_LOCK if called with QF_TRY_LOCK.  */
	uint64_t qf_count_key_value(const QF *qf, uint64_t key, uint64_t value,
															uint8_t flags);

	/* Returns a unique index corresponding to the key in the CQF.  Note
		 that this can change if further modifications are made to the
		 CQF.

		 If the key is not found then returns QF_DOESNT_EXIST.
		 May return QF_COULDNT_LOCK if called with QF_TRY_LOCK.
	 */
	int64_t qf_get_unique_index(const QF *qf, uint64_t key, uint64_t value,
															uint8_t flags);


	/****************************************
   Metadata accessors.
	****************************************/

	/* Hashing info */
	enum qf_hashmode qf_get_hashmode(const QF *qf);
	uint64_t         qf_get_hash_seed(const QF *qf);
	__uint128_t      qf_get_hash_range(const QF *qf);

	/* Space usage info. */
	bool     qf_is_auto_resize_enabled(const QF *qf);
	uint64_t qf_get_total_size_in_bytes(const QF *qf);
	uint64_t qf_get_nslots(const QF *qf);
	uint64_t qf_get_num_occupied_slots(const QF *qf);

	/* Bit-sizes info. */
	uint64_t qf_get_num_key_bits(const QF *qf);
	uint64_t qf_get_num_value_bits(const QF *qf);
	uint64_t qf_get_num_key_remainder_bits(const QF *qf);
	uint64_t qf_get_bits_per_slot(const QF *qf);

	/* Number of (distinct) key-value pairs. */
	uint64_t qf_get_sum_of_counts(const QF *qf);
	uint64_t qf_get_num_distinct_key_value_pairs(const QF *qf);
	
	/****************************************
		Iterators
	*****************************************/
	
	typedef struct quotient_filter_iterator quotient_filter_iterator;
	typedef quotient_filter_iterator QFi;

#define QF_INVALID (-4)
#define QFI_INVALID (-5)
	
	/* Initialize an iterator starting at the given position.
	 * Return value:
	 *  >= 0: iterator is initialized and positioned at the returned slot.
	 *   = QFI_INVALID: iterator has reached end.
	 */
	int64_t qf_iterator_from_position(const QF *qf, QFi *qfi, uint64_t position);

	/* Initialize an iterator and position it at the smallest index
	 * containing a key-value pair whose hash is greater than or equal
	 * to the specified key-value pair.
	 * Return value:
	 *  >= 0: iterator is initialized and position at the returned slot.
	 *   = QFI_INVALID: iterator has reached end.
	 */
	int64_t qf_iterator_from_key_value(const QF *qf, QFi *qfi, uint64_t key,
																		 uint64_t value, uint8_t flags);

	/* Requires that the hash mode of the CQF is INVERTIBLE or NONE.
	 * If the hash mode is DEFAULT then returns QF_INVALID.
	 * Return value:
	 *   = 0: Iterator is still valid.
	 *   = QFI_INVALID: iterator has reached end.
	 *   = QF_INVALID: hash mode is QF_DEFAULT_HASH
	 */
	int qfi_get_key(const QFi *qfi, uint64_t *key, uint64_t *value, uint64_t
									*count);

	/* Return value:
	 *   = 0: Iterator is still valid.
	 *   = QFI_INVALID: iterator has reached end.
	 */
	int qfi_get_hash(const QFi *qfi, uint64_t *hash, uint64_t *value, uint64_t
									 *count);

	/* Advance to next entry.
	 * Return value:
	 *   = 0: Iterator is still valid.
	 *   = QFI_INVALID: iterator has reached end.
	 */
	int qfi_next(QFi *qfi);

	/* Check to see if the if the end of the QF */
	bool qfi_end(const QFi *qfi);

	/************************************
   Miscellaneous convenience functions.
	*************************************/
	
	/* Reset the CQF to an empty filter. */
	void qf_reset(QF *qf);

	/* The caller should call qf_init on the dest QF using the same
	 * parameters as the src QF before calling this function. Note: src
	 * and dest must be exactly the same, including number of slots.  */
	void qf_copy(QF *dest, const QF *src);

	/* merge two QFs into the third one. Note: merges with any existing
		 values in qfc.  */
	void qf_merge(const QF *qfa, const QF *qfb, QF *qfc);

	/* merge multiple QFs into the final QF one. */
	void qf_multi_merge(const QF *qf_arr[], int nqf, QF *qfr);

	/* find cosine similarity between two QFs. */
	uint64_t qf_inner_product(const QF *qfa, const QF *qfb);

	/* square of the L_2 norm of a QF (i.e. sum of squares of counts of
		 all items in the CQF). */
	uint64_t qf_magnitude(const QF *qf);

	/***********************************
		Debugging functions.
	************************************/

	void qf_dump(const QF *);
	void qf_dump_metadata(const QF *qf);


#ifdef __cplusplus
}
#endif

#endif /* _GQF_H_ */
