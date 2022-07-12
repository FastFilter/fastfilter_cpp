
#ifndef TC_SHORTCUT_HPP
#define TC_SHORTCUT_HPP


#include "hashutil.h"
#include "./tc-sym.hpp"

template <typename HashFamily = hashing::TwoIndependentMultiplyShift>
class TC_shortcut {
    HashFamily h0;
    size_t capacity{0};
    const size_t filter_max_capacity;
    //    const size_t remainder_length = 8;
    const size_t quotient_range = tc_sym::QUOTS;
    //    const uint16_t qr_range = 12800UL;
    // const size_t quotient_length = 6;
    const size_t single_pd_capacity = tc_sym::MAX_CAP;

    const size_t number_of_pd;
    const double load;
    __m512i *pd_array;

    size_t insert_existing_counter = 0;
    size_t add_op_counter = 0;

public:
    static inline auto TC_compute_number_of_PD(size_t max_number_of_elements, size_t single_pd_max_cap, double l1_load) -> size_t {
        double b = single_pd_max_cap * l1_load;
        size_t res = (std::size_t) ceil(max_number_of_elements / ((double) b));
        return (res & 1) ? res + 1 : res;
    }
    TC_shortcut(size_t max_number_of_elements, double level1_load_factor)
        : filter_max_capacity(max_number_of_elements),
          number_of_pd(TC_compute_number_of_PD(max_number_of_elements, single_pd_capacity, level1_load_factor)),
          load(level1_load_factor), h0()
    //   pd_index_length(ceil_log2(compute_number_of_PD(max_number_of_elements, max_capacity, level1_load_factor, 1)))
    {
        //        hashing_test = false;
        assert(quotient_range + single_pd_capacity * 9 <= 512);
        assert(!(number_of_pd & 1));
        assert(single_pd_capacity == tc_sym::MAX_CAP);
        int ok = posix_memalign((void **) &pd_array, 64, 64 * number_of_pd);

        if (ok != 0) {
            std::cout << "Failed!!!" << std::endl;
            assert(false);
            return;
        }


        // std::fill(pd_array, pd_array + number_of_pd, __m512i{(INT64_C(1) << 50) - 1, 0, 0, 0, 0, 0, 0, 0});
        static_assert(UINT64_C(-1) == 0xffff'ffff'ffff'ffff);
        // __m512i init_pd = {0xffff'ffff'ffff'ffff, 0x00000000'0000ffff, 0, 0, 0, 0, 0, 0};
        ;
        std::fill(pd_array, pd_array + number_of_pd, __m512i{-1, 0x00000000'0000ffff, 0, 0, 0, 0, 0, 0});
    }

    virtual ~TC_shortcut() {
        free(pd_array);
        // pd_capacity_vec.clear();
    }

    __attribute__((always_inline)) inline static uint32_t reduce32(uint32_t hash, uint32_t n) {
        // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
        return (uint32_t) (((uint64_t) hash * n) >> 32);
    }

    static inline uint32_t reduce(uint64_t hash, uint32_t n) {
        // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
        return (uint32_t) (((hash & 0xffffffffL) * n) >> 32);
    }


    __attribute__((always_inline)) inline static uint16_t fixed_reduce(u32 hash) {
        // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
        constexpr u32 mod = tc_sym::QUOTS << 8;
        // static_assert(mod == 12800UL);
        return (uint16_t) (((uint32_t) hash * mod) >> 16);
    }


    bool v_alt_cfs(uint32_t bin_index, uint32_t rem) const {
        uint32_t alt_bin_index = get_alt_index_cf_stable(bin_index, rem);
        uint32_t alt_alt_bin_index = get_alt_index_cf_stable(alt_bin_index, rem);
        bool res = (bin_index == alt_alt_bin_index);
        if (!res) {
            // auto offset = Morton_offset(rem);
            int diff1 = bin_index - alt_bin_index;
            int diff2 = alt_alt_bin_index - alt_bin_index;
            std::cout << std::string(80, '=') << std::endl;
            std::cout << "number_of_pd: \t" << number_of_pd << std::endl;
            std::cout << std::string(80, '~') << std::endl;
            // std::cout << "offset: \t" << offset << std::endl;
            std::cout << "bin_index: \t" << bin_index << std::endl;
            std::cout << "alt_bin_index: \t" << alt_bin_index << std::endl;
            std::cout << "alt_alt_bin_index: \t" << alt_alt_bin_index << std::endl;
            std::cout << std::string(80, '~') << std::endl;
            std::cout << "diff1: \t" << diff1 << "\t|\t" << (number_of_pd - diff1) << std::endl;
            std::cout << "diff2: \t" << diff2 << "\t|\t" << (number_of_pd - diff2) << std::endl;
            std::cout << std::string(80, '=') << std::endl;

            //            auto offset1 = Morton_offset(rem);
        }
        assert(res);
        return res;
    }

    inline uint32_t get_alt_index_simp(uint32_t bin_index, uint32_t qr) const {
        return (bin_index + qr) % number_of_pd;
    }
    // inline uint32_t get_alt_index_cf_stable(uint32_t bin_index, uint32_t rem) const {
    inline uint32_t get_alt_index_cf_stable(uint32_t bin_index, uint32_t qr) const {
        //Taken from https://github.com/FastFilter/fastfilter_cpp/blob/master/src/cuckoo/cuckoofilter_stable.h
        // (The variables names where changed).
        uint64_t hash = qr * 0xc4ceb9fe1a85ec53L;
        // we don't use xor; instead, we ensure bucketCount is even,
        // and bucket2 = bucketCount - bucket - y,
        // and if negative add the bucketCount,
        // where y is 1..bucketCount - 1 and odd -
        // that way, bucket2 is never the original bucket,
        // and running this twice will give the original bucket, as needed
        uint32_t r = (reduce(hash, number_of_pd >> 1) << 1) + 1;

        // this is needed because the bucket size is not always 2^n:
        int32_t b2 = number_of_pd - bin_index - r;
        if (b2 < 0) {
            b2 += number_of_pd;
        }
        // I tried the following alternatives (also combinations),
        // but performance is the same:

        // uint32_t b2 = bucketCount - index - r;
        // b2 += bucketCount * (b2 >> 31);

        // int32_t b2 = bucketCount - index - r;
        // b2 += bucketCount & (b2 >> 31);

        // int32_t b2 = r - index;
        // b2 += bucketCount & (b2 >> 31);

        return b2;
    }

    inline auto lookup(const uint64_t &item) const -> bool {
        const u64 s = h0(item);
        uint32_t out1 = s >> 32u, out2 = s;
        const uint32_t pd_index1 = reduce32(out1, (uint32_t) number_of_pd);
        const uint16_t qr = fixed_reduce((uint16_t) out2);
        const int64_t quot = qr >> 8;
        const uint8_t rem = qr;

        assert(pd_index1 < number_of_pd);
        assert(quot <= (int64_t)tc_sym::QUOTS);
        // assert(validate_get_alt_index_non_2power(pd_index1, rem));
        assert(v_alt_cfs(pd_index1, qr));
        assert(v_alt_cfs(pd_index1, rem));


        return (tc_sym::find(quot, rem, pd_array + pd_index1)) || (tc_sym::find(quot, rem, pd_array + get_alt_index_cf_stable(pd_index1, qr)));
        //        return (tc_sym::find(quot, rem, pd_array + pd_index1)) || (tc_sym::find(quot, rem, pd_array + get_alt_index_simp(pd_index1, qr)));
        //    (tc_sym::find(quot, rem, pd_array + get_alt_index_non_2power(pd_index1, rem)));
    }

    bool insert_no_shortcut(const uint64_t &item) {
        //        add_op_counter++;
        const u64 s = h0(item);
        uint32_t out1 = s >> 32u, out2 = s;
        const uint32_t pd_index1 = reduce32(out1, (uint32_t) number_of_pd);
        const uint16_t qr = fixed_reduce((uint16_t) out2);
        const int64_t quot = qr >> 8;
        const uint8_t rem = qr;
        assert(pd_index1 < number_of_pd);
        assert(quot <= quotient_range);
        assert(rem <= 255);
        // assert(validate_get_alt_index_non_2power(pd_index1, rem));
        assert(v_alt_cfs(pd_index1, rem));
        assert(v_alt_cfs(pd_index1, qr));

        // const uint32_t pd_index2 = get_alt_index_non_2power(pd_index1, rem);
        const uint32_t pd_index2 = get_alt_index_cf_stable(pd_index1, qr);
        assert(pd_index2 < number_of_pd);

        auto cap1 = tc_sym::get_cap(pd_array + pd_index1);
        auto cap2 = tc_sym::get_cap(pd_array + pd_index2);
        if (cap1 < cap2) {
            auto res = tc_sym::add(quot, rem, &pd_array[pd_index1]);
            assert((!res) or lookup(item));
            return res;

        } else {
            if (cap2 == tc_sym::MAX_CAP) {
                return false;
            }
            auto res = tc_sym::add(quot, rem, &pd_array[pd_index2]);
            assert((!res) or lookup(item));
            return res;
        }
    }

    bool insert(const uint64_t &item) {
        const u64 s = h0(item);
        uint32_t out1 = s >> 32u, out2 = s;
        const uint32_t pd_index1 = reduce32(out1, (uint32_t) number_of_pd);
        const uint16_t qr = fixed_reduce((uint16_t) out2);
        const int64_t quot = qr >> 8;
        const uint8_t rem = qr;
        assert(pd_index1 < number_of_pd);
        assert(quot <= (int64_t)quotient_range);
        assert(rem <= 255);
        // assert(validate_get_alt_index_non_2power(pd_index1, rem));
        assert(v_alt_cfs(pd_index1, rem));
        assert(v_alt_cfs(pd_index1, qr));

        if (tc_sym::pd_less_than_thres(pd_array + pd_index1)) {
#ifndef NDEBUG
            auto res = tc_sym::add(quot, rem, &pd_array[pd_index1]);
            assert((!res) or lookup(item));
            return res;
#endif//!NDEBUG
            return tc_sym::add(quot, rem, &pd_array[pd_index1]);
        }
        const uint32_t pd_index2 = get_alt_index_cf_stable(pd_index1, qr);
        //        const uint32_t pd_index2 = get_alt_index_simp(pd_index1, qr);
        assert(pd_index2 < number_of_pd);

        const uint64_t *word1 = ((const uint64_t *) (pd_array + pd_index1)) + 1;
        const uint64_t *word2 = ((const uint64_t *) (pd_array + pd_index2)) + 1;

        if (word1[0] < word2[0]) {
            assert(tc_sym::get_cap(pd_array + pd_index1) <= tc_sym::get_cap(pd_array + pd_index2));
            if (word1[0] >> 63) {
                assert(tc_sym::pd_full(pd_array + pd_index1));
                return false;
            }
            auto res = tc_sym::add(quot, rem, &pd_array[pd_index1]);
            assert((!res) or lookup(item));
            return res;
        } else {
            assert(tc_sym::get_cap(pd_array + pd_index1) >= tc_sym::get_cap(pd_array + pd_index2));
            if (word2[0] >> 63) {
                assert(tc_sym::pd_full(pd_array + pd_index2));
                return false;
            }
            auto res = tc_sym::add(quot, rem, &pd_array[pd_index2]);
            assert((!res) or lookup(item));
            return res;
        }
    }

    void remove(const uint64_t &item) {
        //        add_op_counter++;
        const u64 s = h0(item);
        assert(lookup(item));
        uint32_t out1 = s >> 32u, out2 = s;
        const uint32_t pd_index1 = reduce32(out1, (uint32_t) number_of_pd);
        const uint16_t qr = fixed_reduce((uint16_t) out2);
        const int64_t quot = qr >> 8;
        const uint8_t rem = qr;
        // assert(validate_get_alt_index_non_2power(pd_index1, rem));
        assert(v_alt_cfs(pd_index1, rem));
        assert(v_alt_cfs(pd_index1, qr));
        tc_sym::conditional_remove(quot, rem, &pd_array[pd_index1]) || tc_sym::remove(quot, rem, &pd_array[get_alt_index_cf_stable(pd_index1, qr)]);
    }

    auto get_capacity() const -> size_t {
        size_t res = 0;
        for (size_t i = 0; i < number_of_pd; ++i) {
            res += tc_sym::get_cap(&pd_array[i]);
        }
        return res;
    }

    size_t get_byte_size() const {
        return number_of_pd * sizeof(pd_array[0]);
    }

    size_t SizeInBytes() const{
        return number_of_pd * sizeof(pd_array[0]);
    }

    auto get_name() const -> std::string {
        return "tc-shortcut";
    }

    auto get_cap() const -> size_t {
        return get_capacity();
    }

    auto get_effective_load() const -> double {
        const size_t slots = number_of_pd * tc_sym::MAX_CAP;
        const size_t filter_cap = get_cap();
        double ratio = 1.0 * filter_cap / slots;
        return ratio;
    }

    std::string info() const {
        std::stringstream ss;
        size_t cap = get_cap();
        ss << std::string(80, '=') << std::endl;
        double ratio = (1.0 * cap) / filter_max_capacity;
        ss << "filter_max_capacity: \t" << filter_max_capacity << std::endl;
        ss << "cap:                 \t" << cap << std::endl;
        ss << "ratio:               \t" << ratio << std::endl;
        ss << "number_of_pd:        \t" << number_of_pd << std::endl;
        ss << "given-load:          \t" << load << std::endl;
        //        ss << "actual-load:         \t" << load << std::endl;
        ss << std::string(80, '=') << std::endl;
        return ss.str();
    }
};

#endif// TWOCHOICER_HEADER
