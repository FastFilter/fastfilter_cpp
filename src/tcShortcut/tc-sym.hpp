/* 
 * My implementation of twoChoicer. A filter with "Power-of-2-choices" insertion policy. 
 * Also known as Vector-Quotient-Filter. 
 * In this version, deletion should work, but was not extensively tested.
 */

#ifndef TC_SYM_HPP
#define TC_SYM_HPP

// #include "../Prefix-Filter/Shift_op.hpp"

#include <algorithm>
#include <assert.h>
#include <immintrin.h>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <stdint.h>
#include <string.h>


typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t u8;

typedef uint_fast64_t uf64;
typedef uint_fast32_t uf32;
typedef uint_fast16_t uf16;
typedef uint_fast8_t uf8;


namespace tc_sym::check {
    template<typename T>
    void zero_array(T *a, size_t a_size) {
        for (size_t i = 0; i < a_size; i++) {
            a[i] = 0;
        }
    }

    auto space_string(std::string s) -> std::string;
    auto to_bin(uint64_t x, size_t length) -> std::string;
    void p_format_word(uint64_t x);
    auto format_word_to_string(uint64_t x, size_t length = 64) -> std::string;

    void print_pd(const __m512i *pd);

    auto validate_number_of_quotient(const __m512i *pd) -> bool;

    auto validate_number_of_quotient(const __m512i *pd, const __m512i *backup_pd) -> bool;
}// namespace tc_sym::check

namespace tc_sym {

    constexpr size_t QUOTS = 80;
    constexpr size_t MAX_CAP = 48;
    constexpr unsigned kBytes2copy = (QUOTS + MAX_CAP + CHAR_BIT - 1) / CHAR_BIT;


    inline uint_fast8_t popcount128(unsigned __int128 x) {
        const uint64_t hi = x >> 64;
        const uint64_t lo = x;
        return _mm_popcnt_u64(lo) + _mm_popcnt_u64(hi);
    }

    inline uint_fast8_t popcount128_on_pd(const __m512i *pd) {
        const auto temp = _mm512_castsi512_si128(*pd);
        return _mm_popcnt_u64(_mm_extract_epi64(temp, 0)) + _mm_popcnt_u64(_mm_extract_epi64(temp, 1));
    }

    __attribute__((always_inline)) inline uint64_t tc_select64(uint64_t x, int64_t j) {
        assert(j < 64);
        const uint64_t y = _pdep_u64(UINT64_C(1) << j, x);
        return _tzcnt_u64(y);
    }

    inline uint64_t select128(unsigned __int128 x, int64_t j) {
        const int64_t pop = _mm_popcnt_u64(x);
        if (j < pop)
            return tc_select64(x, j);
        return 64 + tc_select64(x >> 64, j - pop);
    }

    __attribute__((always_inline)) inline size_t select_on_first_part(int64_t quot, const __m512i *pd) {
        return select128(((unsigned __int128 const *) pd)[0], quot);
    }

    __attribute__((always_inline)) inline size_t select_on_pd_cap(const __m512i *pd) {
        constexpr size_t t = QUOTS - 1;
        u64 Header[2];
        memcpy(Header, pd, 16);
        //  = (const u64*) pd;
        const int64_t pop = _mm_popcnt_u64(Header[0]);
        if ((int64_t)t < pop)
            return tc_select64(Header[0], t);
        return 64 + tc_select64(Header[1], t - pop);
        // return select128(((unsigned __int128 const *) pd)[0], quot);
    }

    __attribute__((always_inline)) inline void select_both_on_word_arr(u64 x, size_t j, uf8 res[2]) {
        assert(j < 64);
        const uint64_t y = _pdep_u64(UINT64_C(3) << j, x);
        assert(_mm_popcnt_u64(y) == 2);

        res[0] = _tzcnt_u64(y);
        res[1] = _tzcnt_u64(_blsr_u64(y));
    }

    inline void select_both_case0_arr(size_t k, u8 res[2], const __m512i *pd) {
        const u64 h0 = ((const u64 *) pd)[0];
        const u64 h1 = ((const u64 *) pd)[1];
        if (k == 0) {
            res[0] = 0;
            res[1] = _tzcnt_u64(h0);
            return;
        }
        auto pop0 = _mm_popcnt_u64(h0);
        if ((int64_t)k < pop0) {
            select_both_on_word_arr(h0, k - 1, res);
            return;
        }
        if ((int64_t)k > pop0) {
            auto new_k = k - pop0;
            select_both_on_word_arr(h1, new_k - 1, res);
            res[0] += 64;
            res[1] += 64;
            return;
        } else {
            res[0] = 63 - _lzcnt_u64(h0);
            res[1] = 64 + _tzcnt_u64(h1);
        }
    }

    /**
     * @brief Get the select mask object
     * 
     * keeps the (j)'th 1 and the (j + 1)'th 1 in x turn on, and turn off the other bits.  
     * @param x 
     * @param j >= 0
     * @return uint64_t 
     */
    inline uint64_t get_select_mask(uint64_t x, int64_t j) {
        assert(_mm_popcnt_u64(x) > j);
        return _pdep_u64(3ul << (j), x);
    }

    /**
     * @brief 
     * 
     * @param x a number with only two set bits. (62 zeros, and 2 ones).
     * @return uint64_t turn on the bits between the currently set bits. 
     *         Also turn off the previously turned on bits. 
     */
    inline uint64_t mask_between_bits(uint64_t x) {
        assert(_mm_popcnt_u64(x) == 2);
        uint64_t hi_bit = (x - 1) & x;
        uint64_t clear_hi = hi_bit - 1;
        uint64_t lo_set = (x - 1);
        uint64_t res = (clear_hi ^ lo_set) & (~x);
        return res;
    }

    inline bool pd_full(const __m512i *pd) {
        uint8_t header_end;
        memcpy(&header_end, reinterpret_cast<const uint8_t *>(pd) + 15, 1);
        return header_end & 128;
    }


    inline uf8 get_cap_naive(const __m512i *pd) {
        auto res = select_on_pd_cap(pd) - (QUOTS - 1);
#ifndef NDEBUG
        auto res0 = select_on_first_part(QUOTS - 1, pd) - (QUOTS - 1);
        assert(res == res0);
#endif//!NDEBUG
        return res;
    }

    __attribute__((always_inline)) inline size_t get_cap(const __m512i *pd) {
        u64 h_last;
        // memcpy(&h_last, (const u8 *) pd + (kBytes2copy - 8), 8);
        memcpy(&h_last, (const u8 *) pd + 8, 8);
        auto temp = _lzcnt_u64(h_last);
        auto res = MAX_CAP - temp;

#ifndef NDEBUG
        auto val = get_cap_naive(pd);
        if (res == val)
            return res;
        std::cout << "res: " << res << std::endl;
        std::cout << "val: " << val << std::endl;

        assert(0);
#endif//!NDEBUG
        assert(res == get_cap_naive(pd));
        assert((res == MAX_CAP) == pd_full(pd));
        return res;
    }

    inline bool pd_less_than_thres(const __m512i *pd) {
        constexpr size_t thres_cap = 36;
        //We need to to test if the last 12+1 bits contains only zeros.
        constexpr u64 thres = (1ULL << (128 - 13 - 64)) - 1;
        const uint64_t h1 = _mm_extract_epi64(_mm512_castsi512_si128(*pd), 1);
        const bool res = h1 < thres;
#ifndef NDEBUG
        auto cap = get_cap(pd);
        bool val = cap < thres_cap;
        if (val != res) {
            get_cap(pd);
            std::cout << "cap: " << cap << std::endl;
            std::cout << "val: " << val << std::endl;
            std::cout << "res: " << res << std::endl;
        }
#endif//!NDEBUG
        return res;
    }

    inline size_t get_spec_quot_cap_inside_word(size_t quot, u64 word) {
        const u64 p_mask = 3 << (quot - 1);
        u64 pdep_res = _pdep_u64(p_mask, word);
        assert(__builtin_popcountll(pdep_res) == 2);
        size_t begin = _tzcnt_u64(pdep_res) + 1;
        size_t end = _tzcnt_u64(_blsr_u64(pdep_res));
        size_t res = end - begin;
        return res;
    }

    inline size_t get_spec_quot_cap(size_t quot, const __m512i *pd) {
        assert(0);
        assert(quot < QUOTS);
        const u64 *pd64 = (const u64 *) pd;
        //        const u64 clean_h = get_clean_header(pd);
        if (quot == 0) {
            std::cout << "h0" << std::endl;
            return _tzcnt_u64(pd64[0]);
        }
        const size_t pop0 = _mm_popcnt_u64(pd64[0]);
        if (quot < pop0) {
            std::cout << "h1" << std::endl;
            return get_spec_quot_cap_inside_word(quot - 1, pd64[0]);
        } else if (pop0 < quot) {
            std::cout << "h2" << std::endl;
            size_t new_q = quot - pop0;
            return get_spec_quot_cap_inside_word(new_q - 1, pd64[1]);
        }
        std::cout << "h3" << std::endl;
        const size_t part1 = _lzcnt_u64(pd64[0]);
        const size_t part2 = _tzcnt_u64(pd64[1]);
        return part1 + part2;
    }

    inline size_t get_spec_quot_cap2(size_t quot, const __m512i *pd) {
        assert(quot < QUOTS);
        const u64 *pd64 = (const u64 *) pd;
        if (quot == 0) {
            //            std::cout << "h0" << std::endl;
            return _tzcnt_u64(pd64[0]);
        }
        size_t begin = select_on_first_part(quot - 1, pd) + 1;
        size_t end = select_on_first_part(quot, pd);
        return end - begin;
    }


    /**
     * @brief Shift 'a' 8 bits right.
     * @param a 
     * @return __m512i 
     */
    inline __m512i shift_right(__m512i a) {
        // The indices encoded here actually encode shift left.
        // idx is an "uint8_t arr[64]" where "arr[i] = i-1" (arr[0] = 0).
        constexpr __m512i idx = {433757350076153919, 1012478732780767239, 1591200115485380623, 2169921498189994007,
                                 2748642880894607391, 3327364263599220775, 3906085646303834159, 4484807029008447543};
        constexpr uint64_t mask = 18446744073709551614ULL;
        __m512i res = _mm512_maskz_permutexvar_epi8(mask, idx, a);
        return res;
    }

    inline __m512i shift_left(__m512i a) {
        constexpr __m512i idx = {578437695752307201, 1157159078456920585, 1735880461161533969, 2314601843866147353, 2893323226570760737, 3472044609275374121, 4050765991979987505, 17801356257212985};

        constexpr uint64_t mask = 9223372036854775807ULL;
        __m512i res = _mm512_maskz_permutexvar_epi8(mask, idx, a);
        return res;
    }

    inline bool pd_find_naive(int64_t quot, uint8_t rem, const __m512i *pd) {
        assert(0 == (reinterpret_cast<uintptr_t>(pd) % 64));
        assert(quot < (int64_t)QUOTS);
        const __m512i target = _mm512_set1_epi8(rem);
        uint64_t v = _mm512_cmpeq_epu8_mask(target, *pd) >> 16ul;

        if (!v) return false;

        // const unsigned __int128 *h = (const unsigned __int128 *) pd;
        // constexpr unsigned __int128 kLeftoverMask = (((unsigned __int128) 1) << (50 + 51)) - 1;
        const unsigned __int128 header = ((const unsigned __int128 *) pd)[0];
        assert(popcount128(header) == QUOTS);
        const uint64_t begin = quot ? select128(header, quot - 1) + 1 - quot : 0;
        const uint64_t end = select128(header, quot) - quot;
        assert(begin <= end);
        assert(end <= MAX_CAP);
        return (v & ((UINT64_C(1) << end) - 1)) >> begin;
    }

    inline bool pd_find_50_v18(int64_t quot, uint8_t rem, const __m512i *pd) {
        assert(0 == (reinterpret_cast<uintptr_t>(pd) % 64));
        assert(quot < (int64_t)QUOTS);
        const __m512i target = _mm512_set1_epi8(rem);
        uint64_t v = _mm512_cmpeq_epu8_mask(target, *pd) >> 16ul;

        if (!v) return false;

        const uint64_t h0 = _mm_extract_epi64(_mm512_castsi512_si128(*pd), 0);
        const uint64_t h1 = _mm_extract_epi64(_mm512_castsi512_si128(*pd), 1);
        if (_blsr_u64(v) == 0) {
            // if ((v << quot)) {
            if ((quot < 64) && (v << quot)) {
                // const unsigned __int128 *h = (const unsigned __int128 *) pd;
                // const unsigned __int128 header = (*h);
                const int64_t mask = v << quot;
                const bool att = (!(h0 & mask)) && (_mm_popcnt_u64(h0 & (mask - 1)) == quot);
                /* if (att != pd_find_naive(quot, rem, pd)) {
                    bool val = pd_find_naive(quot, rem, pd);
                    std::cout << std::string(80, '=') << std::endl;
                    std::cout << "val:  \t" << val << std::endl;
                    std::cout << "att:  \t" << att << std::endl;
                    std::cout << "quot: \t" << quot << std::endl;
                    std::cout << "rem:  \t" << (u16) rem << std::endl;
                    std::cout << "cap:  \t" << get_cap(pd) << std::endl;
                    std::cout << std::string(80, '~') << std::endl;
                    bool a = (!(h0 & mask));
                    bool b = (_mm_popcnt_u64(h0 & (mask - 1)) == quot);
                    std::cout << "a: " << a << std::endl;
                    std::cout << "b: " << b << std::endl;
                    std::cout << std::string(80, '=') << std::endl;

                    assert(att == pd_find_naive(quot, rem, pd));
                }
  */
                assert(att == pd_find_naive(quot, rem, pd));
                return (!(h0 & mask)) && (_mm_popcnt_u64(h0 & (mask - 1)) == quot);
            } else {
                // auto *h = (const unsigned __int128 *) pd;
                // constexpr unsigned __int128 kLeftoverMask = (((unsigned __int128) 1) << (50 + 51)) - 1;
                // const unsigned __int128 header = (*h) & kLeftoverMask;
                const unsigned __int128 header = ((const unsigned __int128 *) pd)[0];
                const unsigned __int128 mask = ((unsigned __int128) v) << quot;
                const bool att = (!(header & mask)) && (popcount128(header & (mask - 1)) == quot);
                assert(att == pd_find_naive(quot, rem, pd));
                return (!(header & mask)) && (popcount128(header & (mask - 1)) == quot);
            }
        }

        const int64_t pop = _mm_popcnt_u64(h0);

        if (quot == 0) {
            // std::cout << "h0" << std::endl;
            return v & (_blsmsk_u64(h0) >> 1ul);
        } else if (quot < pop) {
            // std::cout << "h1" << std::endl;
            const uint64_t mask = (~_bzhi_u64(-1, quot - 1));
            const uint64_t h_cleared_quot_set_bits = _pdep_u64(mask, h0);
            return (((_blsmsk_u64(h_cleared_quot_set_bits) ^ _blsmsk_u64(_blsr_u64(h_cleared_quot_set_bits))) & (~h0)) >> quot) & v;
        } else if (quot > pop) {
            // std::cout << "h2" << std::endl;

            const uint64_t mask = (~_bzhi_u64(-1, quot - pop - 1));
            const uint64_t h_cleared_quot_set_bits = _pdep_u64(mask, h1);
            return (((_blsmsk_u64(h_cleared_quot_set_bits) ^ _blsmsk_u64(_blsr_u64(h_cleared_quot_set_bits))) & (~h1)) >> (quot - pop)) & (v >> (64 - pop));
        } else {
            // std::cout << "h3" << std::endl;

            const uint64_t helper = _lzcnt_u64(h0);
            const uint64_t temp = (63 - helper) + 1;
            const uint64_t diff = helper + _tzcnt_u64(h1);
            return diff && ((v >> (temp - quot)) & ((UINT64_C(1) << diff) - 1));
        }
    }

    inline bool find(int64_t quot, uint8_t rem, const __m512i *pd) {
#ifndef NDEBUG
        bool val = pd_find_naive(quot, rem, pd);
        bool res = pd_find_50_v18(quot, rem, pd);
        if (val != res) {
            std::cout << "val:  \t" << val << std::endl;
            std::cout << "res:  \t" << res << std::endl;
            std::cout << "quot: \t" << quot << std::endl;
            std::cout << "rem:  \t" << (u16) rem << std::endl;
            std::cout << "cap:  \t" << get_cap(pd) << std::endl;
            assert(0);
        }
#endif//!NDEBUG
        return pd_find_50_v18(quot, rem, pd);
    }

    inline void body_remove_case0_avx(size_t index, __m512i *pd) {
        const __m512i shifted_pd = shift_left(*pd);
        const uint64_t mask = _bzhi_u64(-1, kBytes2copy + index);
        _mm512_store_si512(pd, _mm512_mask_blend_epi8(mask, shifted_pd, *pd));
    }

    inline void header_remove_branch(size_t h_end, __m512i *pd) {
#ifndef NDEBUG
        const __m512i backup = *pd;
#endif// DEBUG                                                   \
        // assert(check::validate_number_of_quotient_case0(pd)); \
        // assert(get_cap_wrap(pd));                             \
        // assert(get_cap_wrap(pd) <= MAX_CAP0);
        constexpr uint64_t J = 64;
        uint64_t *pd64 = ((uint64_t *) pd);
        if (h_end >= 64) {
            const size_t rel_index = h_end - 64u;
            assert(rel_index < J);
            const uint64_t index_mask = _bzhi_u64(-1, rel_index);
            uint64_t lo = pd64[1] & index_mask;
            uint64_t hi = ((pd64[1] & ~index_mask) >> 1u);// & _bzhi_u64(-1, J);
            // uint64_t hi = (pd64[1] >> rel_index) << (rel_index + 1);
            assert(!(lo & hi));
            const uint64_t new_h1 = (lo | hi);
            memcpy(pd64 + 1, &new_h1, 8);
//            assert(check::validate_number_of_quotient(pd, &backup));
//            assert(check::validate_number_of_quotient(pd));
            return;
        }
        uint64_t h0, h1;
        // memcpy(&h0, (const uint64_t *) pd, 8);
        memcpy(&h0, pd64, 8);
        memcpy(&h1, pd64 + 1, 8);
        const uint64_t bit_to_move = pd64[1] << 63u;
        const uint64_t h0_mask = _bzhi_u64(-1, h_end);

        const uint64_t h0_lower = h0 & h0_mask;
        const uint64_t h0_hi = (h0 & ~h0_mask) >> 1u;
        assert(!(h0_lower & h0_hi));
        pd64[0] = h0_lower | h0_hi | bit_to_move;

        pd64[1] >>= 1u;
        // const uint64_t new_h1 = (pd64[1] >> 1u
        // memcpy(pd64 + 1, &new_h1, 5);
//        assert(check::validate_number_of_quotient(pd, &backup));
//        assert(check::validate_number_of_quotient(pd));
    }

    inline bool conditional_remove(int64_t quot, uint8_t rem, __m512i *pd) {
        assert(quot < (int64_t)QUOTS);

        const __m512i target = _mm512_set1_epi8(rem);
        uint64_t v = _mm512_cmpeq_epu8_mask(target, *pd) >> kBytes2copy;
        if (!v)
            return false;


        /* if (_blsr_u64(v) == 0) {
            const uint64_t i = _tzcnt_u64(v);

            const bool find_res = (!(header & (((unsigned __int128) 1) << (quot + i)))) &&
                                  (pd512::popcount128(header & (((unsigned __int128) 1 << (quot + i)) - 1)) == quot);

            if (find_res) {
                assert(rem == ((const uint8_t *) pd)[kBytes2copy + i]);

                const uint64_t shift = i + quot;
                unsigned __int128 new_header = header & ((((unsigned __int128) 1) << shift) - 1);
                new_header |= ((header >> (shift + 1)) << (shift));
                // new_header |= ((header >> shift) << (shift - 1));

                assert(pd512::popcount128(header) == 50);

                assert(pd512::validate_number_of_quotient(pd));
                memcpy(pd, &new_header, kBytes2copy);
                assert(pd512::validate_number_of_quotient(pd));

                memmove(&((uint8_t *) pd)[kBytes2copy + i],
                        &((const uint8_t *) pd)[kBytes2copy + i + 1],
                        sizeof(*pd) - (kBytes2copy + i + 1));
            }
            return find_res;
        }
 */

        uf8 res[2] = {0, 0};
        select_both_case0_arr(quot, res, pd);
        const uf64 begin = res[0] + (quot != 0);
        uf64 end = res[1];
        assert(begin <= end);
        const uint64_t begin_fingerprint = begin - quot;
        const uint64_t end_fingerprint = end - quot;

        assert(begin_fingerprint <= end_fingerprint);
        assert(end_fingerprint <= MAX_CAP);

        uint64_t v1 = v >> begin_fingerprint;
        uint64_t i = _tzcnt_u64(v1) + begin_fingerprint;

        // assert(((1 << end_fingerprint) < v) == (i < end_fingerprint));
        if (i >= end_fingerprint)
            return false;

        assert(pd_find_naive(quot, rem, pd));
        assert(rem == ((const uint8_t *) pd)[kBytes2copy + i]);

        const uint64_t shift = i + quot;

        header_remove_branch(shift, pd);
        body_remove_case0_avx(i, pd);

        return true;
    }

    inline bool remove(int64_t quot, uint8_t rem, __m512i *pd) {
        assert(quot < (int64_t)QUOTS);
        assert(find(quot, rem, pd));

        const __m512i target = _mm512_set1_epi8(rem);
        uint64_t v = _mm512_cmpeq_epu8_mask(target, *pd) >> kBytes2copy;
        assert(v);


        if (_blsr_u64(v) == 0) {
            const uint64_t i = _tzcnt_u64(v);
            size_t h_index = i + quot;

            header_remove_branch(h_index, pd);
            body_remove_case0_avx(i, pd);
            return true;
        }

        uf8 res[2] = {0, 0};
        select_both_case0_arr(quot, res, pd);
        const uf64 begin = res[0] + (quot != 0);
        uf64 end = res[1];
        assert(end == select_on_first_part(quot, pd));
        assert(begin <= end);
        const uint64_t begin_fingerprint = begin - quot;
        const uint64_t end_fingerprint = end - quot;

        const uint64_t b_mask = ~((1ull << (begin_fingerprint)) - 1);
        const uint64_t end_mask = ((1ull << end_fingerprint) - 1);
        const uint64_t mask = b_mask & end_mask;
        assert((begin < end) ? mask : !mask);
        const uint64_t v_masked = v & mask;
        assert(v_masked);

        const uint64_t i = _tzcnt_u64(v_masked);

        header_remove_branch(end - 1, pd);
        body_remove_case0_avx(i, pd);

        return true;
    }

    inline void body_add3(size_t end_fingerprint, uint8_t rem, __m512i *pd) {
        const __m512i shifted_pd = shift_right(*pd);
        ((uint8_t *) pd)[kBytes2copy + end_fingerprint] = rem;
        const uint64_t mask = _bzhi_u64(-1, kBytes2copy + end_fingerprint + 1);
        _mm512_store_si512(pd, _mm512_mask_blend_epi8(mask, shifted_pd, *pd));
    }

    inline void write_header6(uint64_t index, __m512i *pd) {
//        assert(check::validate_number_of_quotient(pd));
        // v_pd512_plus::print_headers(pd);
        // constexpr uint64_t h1_mask = ((1ULL << (101 - 64)) - 1);

        uint64_t *pd64 = (uint64_t *) pd;
        // const uint64_t low_h1 = ((pd64[1] << 1) & h1_const_mask) | (pd64[0] >> 63u);
        const uint64_t low_h1 = (pd64[1] << 1) | (pd64[0] >> 63u);

        const uint64_t h0_mask = _bzhi_u64(-1, index);
        const uint64_t h0_lower = pd64[0] & h0_mask;

        // pd64[0] = h0_lower | ((pd64[0] << 1u) & ~h0_mask);
        pd64[0] = h0_lower | ((pd64[0] & ~h0_mask) << 1u);
        memcpy(pd64 + 1, &low_h1, kBytes2copy - 8);
        // v_pd512_plus::print_headers(pd);

//        assert(check::validate_number_of_quotient(pd));
    }


    inline void header_remove_naive(uint64_t index, __m512i *pd) {
//        assert(check::validate_number_of_quotient(pd));

        const unsigned __int128 *h = (const unsigned __int128 *) pd;
        constexpr unsigned __int128 kLeftoverMask = (((unsigned __int128) 1) << (50 + 51)) - 1;
        const unsigned __int128 header = (*h) & kLeftoverMask;

        // const uint64_t end = select128(header, quot);

        const unsigned __int128 mask = (((unsigned __int128) 1) << index) - 1ull;
        unsigned __int128 high_header = (header & ~mask) >> 1ull;
        unsigned __int128 low_header = header & mask;
        // assert(!(high_header & low_header));
        unsigned __int128 new_header = high_header | low_header;
        memcpy(pd, &new_header, kBytes2copy);
//        assert(check::validate_number_of_quotient(pd));
    }

    inline void header_remove(uint64_t index, __m512i *pd) {
//        assert(tc_sym::check::validate_number_of_quotient(pd));
        // v_pd512_plus::print_headers(pd);

        uint64_t *pd64 = (uint64_t *) pd;
        const uint64_t h1_lsb = pd64[1] << 63u;
        // const uint64_t low_h1 = (pd64[1] & h1_const_mask) >> 1u;
        pd64[1] >>= 1u;
        // const uint64_t low_h1 = pd64[1] & h1_const_mask) >> 1u;
        const uint64_t h0_mask = _bzhi_u64(-1, index);
        const uint64_t h0_lower = pd64[0] & h0_mask;
        const uint64_t h0_higher = ((pd64[0] & ~h0_mask) >> 1u) | h1_lsb;

        pd64[0] = h0_lower | h0_higher;
        // memcpy(pd64 + 1, &low_h1, kBytes2copy - 8);

        // v_pd512_plus::print_headers(pd);
//        assert(tc_sym::check::validate_number_of_quotient(pd));
    }


    inline void add_quot0(uint8_t rem, __m512i *pd) {
//        assert(check::validate_number_of_quotient(pd));
        uint64_t *pd64 = (uint64_t *) pd;
        //TODO: end can always be zero.
        const uint64_t end = _tzcnt_u64(pd64[0]);
        const uint64_t low_h1 = ((pd64[1] << 1)) | (pd64[0] >> 63u);
        memcpy(pd64 + 1, &low_h1, kBytes2copy - 8);
        pd64[0] <<= 1u;
//        assert(check::validate_number_of_quotient(pd));

        body_add3(end, rem, pd);
    }


    /**
     * @brief Inserting a new element when pd is not full. 
     * 
     * @param quot 
     * @param rem 
     * @param pd 
     */
    inline void add_wrap(int64_t quot, uint8_t rem, __m512i *pd) {
        assert(!pd_full(pd));
        assert(quot < (int64_t)QUOTS);

        if (quot == 0) {
            add_quot0(rem, pd);
            return;
        }


        // constexpr uint64_t h1_const_mask = ((1ULL << (101 - 64)) - 1);

        const uint64_t h0 = _mm_extract_epi64(_mm512_castsi512_si128(*pd), 0);
        const uint64_t pop = _mm_popcnt_u64(h0);
        // const uint64_t pop = (h0 ^ (h0 >> 8u)) % QUOTS;
        if (quot < (int64_t)pop) {
            const uint64_t end = tc_select64(h0, quot);
            write_header6(end, pd);
            // write_header_naive(quot, end, pd);
            body_add3(end - quot, rem, pd);
            return;
        } else {
            const uint64_t h1 = _mm_extract_epi64(_mm512_castsi512_si128(*pd), 1);
            const uint64_t end = tc_select64(h1, quot - pop);

            //header
//            assert(tc_sym::check::validate_number_of_quotient(pd));

            const uint64_t h1_mask = _bzhi_u64(-1, end);
            const uint64_t h1_low = h1 & h1_mask;
            assert(end < 64);
            const uint64_t h1_high = (h1 >> end) << (end + 1);
            assert(!(h1_low & h1_high));
            const uint64_t new_h1 = h1_low | h1_high;

            memcpy(&((uint64_t *) pd)[1], &new_h1, kBytes2copy - 8);
//            assert(tc_sym::check::validate_number_of_quotient(pd));

            //body
            const uint64_t index = 64 + end - quot;
            // body_add_naive(index, rem, pd);
            body_add3(index, rem, pd);

//            assert(tc_sym::check::validate_number_of_quotient(pd));
            assert(find(quot, rem, pd));

            return;
        }
    }


    inline bool add(int64_t quot, uint8_t rem, __m512i *pd) {
        // return add_db(quot, rem, pd);
        assert(!pd_full(pd));
        add_wrap(quot, rem, pd);
        assert(find(quot, rem, pd));
        return true;
    }
}// namespace tc_sym


#endif// TC_SYM_HPP
