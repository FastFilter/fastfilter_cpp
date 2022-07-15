#ifndef FILTERS_MIN_PD256_HPP
#define FILTERS_MIN_PD256_HPP

#include <cassert>
#include <climits>

#include <iostream>

#include <immintrin.h>
#include <x86intrin.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t u8;



namespace min_pd {
    constexpr size_t QUOTS = 25;
    constexpr size_t MAX_CAP0 = 25;
    constexpr u64 H_mask = (1ULL << (QUOTS + MAX_CAP0)) - 1;
    struct add_res {
        u8 quot;
        u8 rem;
        bool passed;
    };

    bool find_core(int64_t quot, uint8_t rem, const __m256i *pd);

    /**
     * Remark: if `x` has `j` set bits, then: select(x,j) == 64. 
     * @param x 
     * @param j 
     * @return The position (starting from 0) of the jth set bit of x.
     */
    inline uint64_t pd_select64(uint64_t x, int64_t j) {
        assert(j < 64);
        const uint64_t y = _pdep_u64(UINT64_C(1) << j, x);
        return _tzcnt_u64(y);
    }

    inline unsigned get_status(const __m256i *pd) {
        u8 temp;
        memcpy(&temp, pd, 1);
        return temp & 31;
    }

    inline unsigned get_h_first_byte(const __m256i *pd) {
        u8 temp;
        memcpy(&temp, pd, 1);
        return temp;
    }

    inline uint8_t get_last_byte(const __m256i *pd) {
        uint8_t x;
        memcpy(&x, ((uint8_t *) pd) + 31, 1);
        return x;
    }

    inline uint64_t get_header(const __m256i *pd) {
        return ((uint64_t *) pd)[0] & ((1ULL << 56) - 1);
    }

    inline uint64_t get_clean_header(const __m256i *pd) {
        // uint64_t res = (_mm_cvtsi128_si64(_mm256_castsi256_si128(*pd)) >> 6) & ((1ULL << 50) - 1);
        return ((((uint64_t *) pd)[0]) >> 6ul) & ((1ULL << 50ul) - 1);
    }

    inline bool did_pd_overflowed(const __m256i *pd) {
        return !(((uint64_t *) pd)[0] & 32);
    }

    inline void set_overflow_bit(__m256i *pd) {
        uint64_t *h_array = ((uint64_t *) pd);
        h_array[0] |= 32;
        h_array[0] ^= 32;
        assert(did_pd_overflowed(pd));
    }

    inline void clear_overflow_bit(__m256i *pd) {
        uint64_t *h_array = ((uint64_t *) pd);
        h_array[0] |= 32;
        // h_array[0] ^= 32;
        assert(!did_pd_overflowed(pd));
    }

    inline uint64_t decode_last_quot(const __m256i *pd) {
        return ((uint64_t *) pd)[0] & 31;
    }

    inline size_t get_cap_naive(const __m256i *pd) {
        // constexpr u64 mask = (1ULL << (QUOTS + MAX_CAP0)) - 1;
        constexpr size_t t = QUOTS - 1;
        const uint64_t header = reinterpret_cast<const u64 *>(pd)[0];
        u64 h0 = (header >> 6);// & H_mask;
        size_t res = pd_select64(h0, t) - t;
        return res;
        //        auto pop0 = _mm_popcnt_u64(h0);
        //        assert(pop0 == QUOTS);
        //        size_t zeros = 64 - pop0;
        //        size_t cap = zeros - 14;
        //        return cap;
    }

    inline size_t get_capacity(const __m256i *pd) {
        const uint64_t header = reinterpret_cast<const u64 *>(pd)[0];
        auto temp = _lzcnt_u64(header << 8);
        uint64_t res = MAX_CAP0 - temp;
        assert(res == get_cap_naive(pd));
        return res;
    }

    inline size_t get_cap(const __m256i *pd) {
        const uint64_t header = reinterpret_cast<const u64 *>(pd)[0];
        auto temp = _lzcnt_u64(header << 8);
        uint64_t res = MAX_CAP0 - temp;
        assert(res == get_cap_naive(pd));
        return res;
    }

    inline bool is_pd_full(const __m256i *pd) {
        const uint64_t header = reinterpret_cast<const u64 *>(pd)[0];
        bool res = header & (1ULL << 55);
        assert(res == (get_cap(pd) == MAX_CAP0));
        return res;
    }

    inline bool pd_full(const __m256i *pd) {
        return is_pd_full(pd);
    }

    inline bool is_header_full(const __m256i *pd) {
        u8 temp;
        memcpy(&temp, (const u8 *) pd + 6, 1);
        return temp & 128;
    }

    inline unsigned header_last_two_bits(const __m256i *pd) {
        const uint64_t header = reinterpret_cast<const u64 *>(pd)[0];
        unsigned res = (header >> 54) & 3;
        return res;
    }

    inline __m256i shift_right(__m256i a) {
        constexpr __m256i idx = {0x605040302010000, 0xe0d0c0b0a090807, 0x161514131211100f, 0x1e1d1c1b1a191817};
        // constexpr uint64_t mask = 18446744073709551614ULL;
        __m256i res = _mm256_permutexvar_epi8(idx, a);
        return res;
    }


    inline size_t get_spec_quot_cap(size_t quot, const __m256i *pd) {
        assert(quot < QUOTS);
        const u64 clean_h = get_clean_header(pd);
        if (quot == 0) {
            return _tzcnt_u64(clean_h);
        }
        const u64 p_mask = 3 << (quot - 1);
        u64 pdep_res = _pdep_u64(p_mask, clean_h);
        assert(__builtin_popcountll(pdep_res) == 2);
        size_t begin = _tzcnt_u64(pdep_res) + 1;
        size_t end = _tzcnt_u64(_blsr_u64(pdep_res));
        size_t res = end - begin;
        return res;
    }

    inline void body_add_case0_avx(size_t body_index, uint8_t rem, __m256i *pd) {
        constexpr unsigned kBytes2copy = 7;
        const __m256i shifted_pd = shift_right(*pd);
        ((uint8_t *) pd)[kBytes2copy + body_index] = rem;
        const uint64_t mask = _bzhi_u64(-1, kBytes2copy + body_index + 1);
        _mm256_store_si256(pd, _mm256_mask_blend_epi8(mask, shifted_pd, *pd));
    }
    inline void body_add_simple(size_t body_index, uint8_t rem, __m256i *pd) {
        // const size_t body_index = end - quot;
        auto mp = (u8 *) pd + 7 + body_index;
        const size_t b2m = (32 - 7) - (body_index + 1);
        memmove(mp + 1, mp, b2m);
        mp[0] = rem;
    }

    inline bool body_add(size_t body_index, u8 rem, __m256i *pd) {
        __m256i pd0 = *pd;
        __m256i pd1 = *pd;
        body_add_case0_avx(body_index, rem, &pd0);
        body_add_simple(body_index, rem, &pd1);
        bool res = memcmp(&pd0, &pd1, 32) == 0;
        if (res)
            return true;

        // std::cout << "h0:\t" << str_bitsMani::format_word_to_string(((u64 *) &pd0)[0]) << std::endl;
        // std::cout << "h1:\t" << str_bitsMani::format_word_to_string(((u64 *) &pd1)[0]) << std::endl;
        // std::cout << "h1:\t" << str_bitsMani::format_word_to_string(((u64 *) &pd1)[0]) << std::endl;
        // std::cout << "h1:\t" << str_bitsMani::format_word_to_string(((u64 *) &pd1)[0]) << std::endl;
        assert(0);
        return false;
    }

    inline size_t get_last_occupied_quot_only_full_pd(const __m256i *pd) {
        assert(is_pd_full(pd));
        const u64 clean_h = get_clean_header(pd);
        const size_t last_quot = pd_select64(~clean_h, MAX_CAP0 - 1) - (MAX_CAP0 - 1);
        assert(last_quot < QUOTS);
#ifndef NDEBUG
        for (size_t i = last_quot + 1; i < QUOTS; i++) {
            assert(!get_spec_quot_cap(i, pd));
        }
#endif//!NDEBUG
        return last_quot;
    }


    inline void sort_k_last_rem(size_t k, __m256i *pd) {
        //   taken from https://stackoverflow.com/questions/2786899/fastest-sort-of-fixed-length-6-int-array/2789530#2789530
        int i, j;
        u8 *d = (u8 *) pd + 32 - k;
        for (i = 1; i < (int) k; i++) {
            u8 tmp = d[i];
            for (j = i; j >= 1 && tmp < d[j - 1]; j--)
                d[j] = d[j - 1];
            d[j] = tmp;
        }
    }

    /**
     * @brief Sorts the remainders in the last run in non-decsending order, and returns the value of the last quot. 
     * 
     * @param pd 
     * @return size_t 
     */
    inline size_t sort_last_quot(__m256i *pd) {
        assert(is_pd_full(pd));
        const u64 clean_h = get_clean_header(pd);
        const size_t last_quot = pd_select64(~clean_h, MAX_CAP0 - 1) - (MAX_CAP0 - 1);
        assert(get_last_occupied_quot_only_full_pd(pd) == last_quot);
        size_t last_zero_index = last_quot + (MAX_CAP0 - 1);
        u64 shifted_h = clean_h << (63 - last_zero_index);
        assert(!_bextr_u64(shifted_h, 63, 1));
        const size_t lq_cap = _lzcnt_u64(shifted_h);
        assert(lq_cap == get_spec_quot_cap(last_quot, pd));
        sort_k_last_rem(lq_cap, pd);
        return last_quot;
    }


    inline size_t get_last_occ_quot_cap(const __m256i *pd) {
        assert(is_pd_full(pd));
        const u64 dirty_h = reinterpret_cast<const u64 *>(pd)[0] >> 6;
        const size_t last_quot = pd_select64(~dirty_h, MAX_CAP0 - 1) - (MAX_CAP0 - 1);
        assert(get_last_occupied_quot_only_full_pd(pd) == last_quot);
        size_t last_zero_index = last_quot + (MAX_CAP0 - 1);
        u64 shifted_h = dirty_h << (63 - last_zero_index);
        assert(!_bextr_u64(shifted_h, 63, 1));
        size_t lq_cap = _lzcnt_u64(shifted_h);
        assert(lq_cap == get_spec_quot_cap(last_quot, pd));
        return lq_cap;
    }


    inline void update_status_old(int64_t last_quot, __m256i *pd) {
        //check this
        uint8_t byte_to_write = last_quot | (get_header(pd) & (32 + 64 + 128));
        memcpy(pd, &byte_to_write, 1);
        assert((int64_t)decode_last_quot(pd) == last_quot);
    }

    inline void update_status(size_t last_occ_quot, __m256i *pd) {
        assert(last_occ_quot < 32);
        auto pd64 = reinterpret_cast<u64 *>(pd);
        pd64[0] = ((pd64[0] | 31) ^ 31) | last_occ_quot;
    }


    inline void pd_add_50_only_rem(uint8_t rem, size_t quot_capacity, __m256i *pd) {
        //FIXME use cmp_leq_as mask for shift avx.
        // constexpr unsigned kBytes2copy = 7;

        // const __m256i target = _mm256_set1_epi8(rem);
        const uint64_t begin_fingerprint = MAX_CAP0 - quot_capacity;
        const uint64_t end_fingerprint = MAX_CAP0;
        assert(begin_fingerprint < end_fingerprint);
        assert(end_fingerprint <= MAX_CAP0);

        uint64_t i = begin_fingerprint;
        auto cp0 = (const u8 *) pd + 7;
        for (; i < end_fingerprint; ++i) {
            if (rem <= cp0[i]) break;
        }
        assert((i == end_fingerprint) ||
               (rem <= ((const uint8_t *) pd)[7 + i]));

        u8 *mp = (u8 *) pd + 7 + i;
        size_t b2m = (32 - 7) - (i + 1);
        assert(b2m < 32);
        memmove(mp + 1, mp, b2m);
        mp[0] = rem;
    }


    inline void add_full_pd(bool did_pd_previously_overflowed, size_t last_occ_quot, int64_t quot, u8 rem, __m256i *pd) {
        assert(is_pd_full(pd));
        assert(last_occ_quot < QUOTS);
        assert(get_spec_quot_cap(last_occ_quot, pd));
#ifndef NDEBUG
        size_t v_lq_cap = get_spec_quot_cap(last_occ_quot, pd);
#endif//!NDEBUG
        const uint64_t dirty_h0 = reinterpret_cast<const u64 *>(pd)[0];
        const size_t set_last_zero_index = (MAX_CAP0 + last_occ_quot + 6);// -1?
        const u64 set_last_zero = 1ULL << set_last_zero_index;
        const bool change_last_quot = (!did_pd_previously_overflowed) | (dirty_h0 & (set_last_zero >> 2));
#ifndef NDEBUG
        if (did_pd_previously_overflowed)
            assert(change_last_quot == (get_spec_quot_cap(last_occ_quot, pd) == 1));
#endif//!NDEBUG
        const size_t end = pd_select64(dirty_h0 >> 6, quot);
        // constexpr u64 h_mask = ((1ULL << 56) - 1);

        const size_t h_index = end + 6;
        const u64 mask = _bzhi_u64(-1, h_index);
        const u64 lo = dirty_h0 & mask;
        // const u64 pre_hi = ((dirty_h0 & ~mask) << 1u);            // | set_last_zero;// & h_mask;
        const u64 hi = ((dirty_h0 & ~mask) << 1u) | set_last_zero;// & h_mask;


        assert(!(lo & hi));
        const u64 h7 = lo | hi;
        memcpy(pd, &h7, 7);

        const size_t body_index = end - quot;
        body_add_case0_avx(body_index, rem, pd);

        if (!change_last_quot)
            return;

        size_t new_last_occ_quot = sort_last_quot(pd);
        update_status(new_last_occ_quot, pd);
    }


    inline add_res new_pd_swap_short(int64_t quot, uint8_t rem, __m256i *pd) {
        uint64_t last_quot = decode_last_quot(pd);
        const bool did_ovf = did_pd_overflowed(pd);
        if (!did_ovf) {
            last_quot = sort_last_quot(pd);
            auto pd64 = reinterpret_cast<u64 *>(pd);
            pd64[0] ^= (last_quot | 32);
        }

        if (last_quot < (uint64_t) quot) {
//            assert(check::val_last_quot_is_sorted(pd));
            return {(u8) quot, rem, false};
        } else if ((int64_t)last_quot == quot) {
            const uint8_t old_rem = get_last_byte(pd);
            if (old_rem < rem) {//todo can be <= instead of <.
                return {(u8) quot, rem, false};
            }

            size_t quot_capacity = get_spec_quot_cap(quot, pd);
            assert(get_spec_quot_cap(quot, pd));
            pd_add_50_only_rem(rem, quot_capacity, pd);
            assert(find_core(quot, rem, pd));

            return {(u8) last_quot, old_rem, false};
        } else {
            const u8 old_rem = get_last_byte(pd);
            add_full_pd(did_ovf, last_quot, quot, rem, pd);

            assert(find_core(quot, rem, pd));
            return {(u8) last_quot, old_rem, false};
        }
    }


    inline bool find_core(int64_t quot, uint8_t rem, const __m256i *pd) {
        assert(0 == (reinterpret_cast<uintptr_t>(pd) % 32));
        assert(quot < (int64_t)QUOTS);

        const __m256i target = _mm256_set1_epi8(rem);
        const uint64_t v = _mm256_cmpeq_epu8_mask(target, *pd) >> 7ul;
        if (!v) return false;

        const uint64_t v_off = _blsr_u64(v);
        const uint64_t h0 = get_clean_header(pd);

        if (v_off == 0) {
            const uint64_t mask = v << quot;
            return (_mm_popcnt_u64(h0 & (mask - 1)) == quot) && (!(h0 & mask));
        }

        if (quot == 0)
            return v & (_blsmsk_u64(h0) >> 1ul);

        uint64_t new_v = (v << quot) & ~h0;
        const uint64_t mask = (~_bzhi_u64(-1, quot - 1));
        const uint64_t h_cleared_quot_set_bits = _pdep_u64(mask, h0);
        const uint64_t h_cleared_quot_plus_one_set_bits = _blsr_u64(h_cleared_quot_set_bits);
        const uint64_t v_mask = _blsmsk_u64(h_cleared_quot_set_bits) ^ _blsmsk_u64(h_cleared_quot_plus_one_set_bits);
        // bool att = v_mask & new_v;
        return v_mask & new_v;
    }


    /**
     * @brief indicate if we should solely look for the element in the next level. 
     * 
     * @param qr 
     * @param pd 
     * @return true Look only in the second level.
     * @return false Look only in this PD.
     */
    inline bool cmp_qr1(uint16_t qr, const __m256i *pd) {
        if (((uint64_t *) pd)[0] & 32) {
            assert(!did_pd_overflowed(pd));
            return false;
        }
        const uint64_t hfb = ((uint64_t *) pd)[0] & 31;
        const uint16_t old_qr = (hfb << 8ul) | get_last_byte(pd);
        return old_qr < qr;
    }

}// namespace min_pd


#endif// FILTERS_PD256_PLUS_HPP
