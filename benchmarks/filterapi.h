#ifndef FILTERAPI_H
#define FILTERAPI_H
#include <climits>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <vector>
#include <set>
#include <stdio.h>

// morton
#include "compressed_cuckoo_filter.h"
#include "morton_sample_configs.h"
#include "cuckoofilter.h"
#include "cuckoofilter_stable.h"
#include "cuckoo_fuse.h"
#include "xorfilter.h"
#include "xorfilter_10bit.h"
#include "xorfilter_13bit.h"
#include "xorfilter_10_666bit.h"
#include "xorfilter_2.h"
#include "xorfilter_2n.h"
#include "xorfilter_plus.h"
#include "xorfilter_singleheader.h"
#include "fusefilter_singleheader.h"
#include "binaryfusefilter_singleheader.h"
#include "xor_binary_fuse_filter.h"
#include "bloom.h"
#include "counting_bloom.h"
#include "gcs.h"
#ifdef __AVX2__
#include "gqf_cpp.h"
#include "simd-block.h"
#endif
#include "simd-block-fixed-fpp.h"

using namespace std;
using namespace hashing;
using namespace cuckoofilter;
using namespace cuckoofusefilter;
using namespace xorfilter;
using namespace xorfilter2;
using namespace xorfilter2n;
using namespace xorfilter_plus;
using namespace bloomfilter;
using namespace counting_bloomfilter;
using namespace gcsfilter;
using namespace CompressedCuckoo; // Morton filter namespace
#ifdef __AVX2__
using namespace gqfilter;
#endif



// Inlining the "contains" which are executed within a tight loop can be both
// very detrimental or very beneficial, and which ways it goes depends on the
// compiler. It is unclear whether we want to benchmark the inlining of Contains,
// as it depends very much on how "contains" is used. So it is maybe reasonable
// to benchmark it without inlining.
//
#define CONTAIN_ATTRIBUTES  __attribute__ ((noinline))

template<typename Table>
struct FilterAPI {};

template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
struct FilterAPI<CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>> {
  using Table = CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table * table) {
    if (0 != table->Add(key)) {
      throw logic_error("The filter is too small to hold all of the elements");
    }
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Delete(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
struct FilterAPI<CuckooFilterStable<ItemType, bits_per_item, TableType, HashFamily>> {
  using Table = CuckooFilterStable<ItemType, bits_per_item, TableType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table * table) {
    if (0 != table->Add(key)) {
      throw logic_error("The filter is too small to hold all of the elements");
    }
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Delete(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<CuckooFuseFilter<ItemType, FingerprintType>> {
  using Table = CuckooFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table * table) {
    if (0 != table->Add(key)) {
      throw logic_error("The filter is too small to hold all of the elements");
    }
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Delete(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

#ifdef __aarch64__
template <typename HashFamily>
struct FilterAPI<SimdBlockFilterFixed<HashFamily>> {
  using Table = SimdBlockFilterFixed<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(ceil(add_count * 8.0 / CHAR_BIT));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};

#endif

#ifdef __AVX2__
template <typename HashFamily>
struct FilterAPI<SimdBlockFilter<HashFamily>> {
  using Table = SimdBlockFilter<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(ceil(log2(add_count * 8.0 / CHAR_BIT)));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};

template <typename HashFamily>
struct FilterAPI<SimdBlockFilterFixed64<HashFamily>> {
  using Table = SimdBlockFilterFixed64<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(ceil(add_count * 8.0 / CHAR_BIT));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};

template <typename HashFamily>
struct FilterAPI<SimdBlockFilterFixed<HashFamily>> {
  using Table = SimdBlockFilterFixed<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(ceil(add_count * 8.0 / CHAR_BIT));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};
#endif

#ifdef __SSE41__
template <typename HashFamily>
struct FilterAPI<SimdBlockFilterFixed16<HashFamily>> {
  using Table = SimdBlockFilterFixed16<HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(ceil(add_count * 8.0 / CHAR_BIT));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};
#endif

template <typename ItemType, typename FingerprintType>
struct FilterAPI<XorFilter<ItemType, FingerprintType>> {
  using Table = XorFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_naive::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_naive::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};


template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_sorted::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_sorted::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_partiallysorted::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_partiallysorted::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_onehash::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_onehash::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_lowmem::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_lowmem::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_fixedsorted::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_fixedsorted::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }

};
template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_naive4wise::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_naive4wise::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};
template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_partiallysorted4wise::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_partiallysorted4wise::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};
template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_lowmem4wise::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_lowmem4wise::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};


template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_4wise_prefetched::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_4wise_prefetched::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType>
struct FilterAPI<xorbinaryfusefilter_prefetched::XorBinaryFuseFilter<ItemType, FingerprintType>> {
  using Table = xorbinaryfusefilter_prefetched::XorBinaryFuseFilter<ItemType, FingerprintType>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};



class MortonFilter {
    Morton3_8* filter;
    size_t size;
public:
    MortonFilter(const size_t size) {
        filter = new Morton3_8((size_t) (size / 0.95) + 64);
        // filter = new Morton3_8((size_t) (2.1 * size) + 64);
        this->size = size;
    }
    ~MortonFilter() {
        delete filter;
    }
    void Add(uint64_t key) {
        filter->insert(key);
    }
    void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end) {
        size_t size = end - start;
        ::std::vector<uint64_t> k(size);
        ::std::vector<bool> status(size);
        for (size_t i = start; i < end; i++) {
            k[i - start] = keys[i];
        }
        // TODO return value and status is ignored currently
        filter->insert_many(k, status, size);
    }
    inline bool Contain(uint64_t &item) {
        return filter->likely_contains(item);
    };
    size_t SizeInBytes() const {
        // according to morton_sample_configs.h:
        // Morton3_8 - 3-slot buckets with 8-bit fingerprints: 11.7 bits/item
        // (load factor = 0.95)
        // so in theory we could just hardcode the size here,
        // and don't measure it
        // return (size_t)((size * 11.7) / 8);

        return filter->SizeInBytes();
    }
};

template<>
struct FilterAPI<MortonFilter> {
    using Table = MortonFilter;
    static Table ConstructFromAddCount(size_t add_count) {
        return Table(add_count);
    }
    static void Add(uint64_t key, Table* table) {
        table->Add(key);
    }
    static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
        table->AddAll(keys, start, end);
    }
    static void Remove(uint64_t, Table *) {
        throw std::runtime_error("Unsupported");
    }
    CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, Table * table) {
        return table->Contain(key);
    }
};

class XorSingle {
public:
    xor8_s filter; // let us expose the struct. to avoid indirection
    explicit XorSingle(const size_t size) {
        if (!xor8_allocate(size, &filter)) {
            throw ::std::runtime_error("Allocation failed");
        }
    }
    ~XorSingle() {
        xor8_free(&filter);
    }
    bool AddAll(const uint64_t* data, const size_t start, const size_t end) {
        return xor8_buffered_populate(data + start, end - start, &filter);
    }
    inline bool Contain(uint64_t &item) const {
        return xor8_contain(item, &filter);
    }
    inline size_t SizeInBytes() const {
        return xor8_size_in_bytes(&filter);
    }
    XorSingle(XorSingle && o) : filter(o.filter)  {
        o.filter.fingerprints = nullptr; // we take ownership for the data
    }
private:
    XorSingle(const XorSingle & o) = delete;
};


template<>
struct FilterAPI<XorSingle> {
    using Table = XorSingle;
    static Table ConstructFromAddCount(size_t add_count) {
        return Table(add_count);
    }
    static void Add(uint64_t, Table*) {
        throw std::runtime_error("Unsupported");
    }
    static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
        table->AddAll(keys.data(), start, end);
    }
    static void Remove(uint64_t, Table *) {
        throw std::runtime_error("Unsupported");
    }
    CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
        // some compilers are not smart enough to do the inlining properly
        return xor8_contain(key, & table->filter);
    }
};


class FuseSingle {
public:
    fuse8_t filter; // let us expose the struct. to avoid indirection
    explicit FuseSingle(const size_t size) {
        if (!fuse8_allocate(size, &filter)) {
            throw ::std::runtime_error("Allocation failed");
        }
    }
    ~FuseSingle() {
        fuse8_free(&filter);
    }
    bool AddAll(const uint64_t* data, const size_t start, const size_t end) {
        return fuse8_populate(data + start, end - start, &filter);
    }
    inline bool Contain(uint64_t &item) const {
        return fuse8_contain(item, &filter);
    }
    inline size_t SizeInBytes() const {
        return fuse8_size_in_bytes(&filter);
    }
    FuseSingle(FuseSingle && o) : filter(o.filter)  {
        o.filter.fingerprints = nullptr; // we take ownership for the data
    }
private:
    FuseSingle(const FuseSingle & o) = delete;
};


template<>
struct FilterAPI<FuseSingle> {
    using Table = FuseSingle;
    static Table ConstructFromAddCount(size_t add_count) {
        return Table(add_count);
    }
    static void Add(uint64_t, Table*) {
        throw std::runtime_error("Unsupported");
    }
    static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
        table->AddAll(keys.data(), start, end);
    }
    static void Remove(uint64_t, Table *) {
        throw std::runtime_error("Unsupported");
    }
    CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
        // some compilers are not smart enough to do the inlining properly
        return fuse8_contain(key, & table->filter);
    }
};

class BinaryFuseSingle {
public:
    binary_fuse8_t filter; // let us expose the struct. to avoid indirection
    explicit BinaryFuseSingle(const size_t size) {
        if (!binary_fuse8_allocate(size, &filter)) {
            throw ::std::runtime_error("Allocation failed");
        }
    }
    ~BinaryFuseSingle() {
        binary_fuse8_free(&filter);
    }
    bool AddAll(const uint64_t* data, const size_t start, const size_t end) {
        return binary_fuse8_populate(data + start, end - start, &filter);
    }
    inline bool Contain(uint64_t &item) const {
        return binary_fuse8_contain(item, &filter);
    }
    inline size_t SizeInBytes() const {
        return binary_fuse8_size_in_bytes(&filter);
    }
    BinaryFuseSingle(BinaryFuseSingle && o) : filter(o.filter)  {
        o.filter.Fingerprints = nullptr; // we take ownership for the data
    }
private:
    BinaryFuseSingle(const BinaryFuseSingle & o) = delete;
};


template<>
struct FilterAPI<BinaryFuseSingle> {
    using Table = BinaryFuseSingle;
    static Table ConstructFromAddCount(size_t add_count) {
        return Table(add_count);
    }
    static void Add(uint64_t , Table* ) {
        throw std::runtime_error("Unsupported");
    }
    static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
        table->AddAll(keys.data(), start, end);
    }
    static void Remove(uint64_t , Table * ) {
        throw std::runtime_error("Unsupported");
    }
    CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
        // some compilers are not smart enough to do the inlining properly
        return binary_fuse8_contain(key, & table->filter);
    }
};


template<size_t blocksize, int k, typename HashFamily>
struct FilterAPI<SimpleBlockFilter<blocksize,k,HashFamily>> {
  using Table = SimpleBlockFilter<blocksize,k,HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) {
    Table ans(ceil(add_count * 8.0 / CHAR_BIT));
    return ans;
  }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Find(key);
  }
};


template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<XorFilter<ItemType, FingerprintType, HashFamily>> {
  using Table = XorFilter<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};


template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<naive::XorFilter<ItemType, FingerprintType, HashFamily>> {
  using Table = naive::XorFilter<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<prefetch::XorFilter<ItemType, FingerprintType, HashFamily>> {
  using Table = prefetch::XorFilter<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename FingerprintStorageType, typename HashFamily>
struct FilterAPI<XorFilter2<ItemType, FingerprintType, FingerprintStorageType, HashFamily>> {
  using Table = XorFilter2<ItemType, FingerprintType, FingerprintStorageType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename HashFamily>
struct FilterAPI<XorFilter10<ItemType, HashFamily>> {
  using Table = XorFilter10<ItemType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename HashFamily>
struct FilterAPI<XorFilter13<ItemType, HashFamily>> {
  using Table = XorFilter13<ItemType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename HashFamily>
struct FilterAPI<XorFilter10_666<ItemType, HashFamily>> {
  using Table = XorFilter10_666<ItemType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename FingerprintStorageType, typename HashFamily>
struct FilterAPI<XorFilter2n<ItemType, FingerprintType, FingerprintStorageType, HashFamily>> {
  using Table = XorFilter2n<ItemType, FingerprintType, FingerprintStorageType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, typename FingerprintType, typename HashFamily>
struct FilterAPI<XorFilterPlus<ItemType, FingerprintType, HashFamily>> {
  using Table = XorFilterPlus<ItemType, FingerprintType, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<GcsFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = GcsFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t, Table*) {
    throw std::runtime_error("Unsupported");
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

#ifdef __AVX2__
template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<GQFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = GQFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};
#endif

template <typename ItemType, size_t bits_per_item, bool branchless, typename HashFamily>
struct FilterAPI<BloomFilter<ItemType, bits_per_item, branchless, HashFamily>> {
  using Table = BloomFilter<ItemType, bits_per_item, branchless, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys.data(), start, end);
  }
  static void Remove(uint64_t, Table *) {
    throw std::runtime_error("Unsupported");
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, bool branchless, typename HashFamily>
struct FilterAPI<CountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily>> {
  using Table = CountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, bool branchless, typename HashFamily>
struct FilterAPI<SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily>> {
  using Table = SuccinctCountingBloomFilter<ItemType, bits_per_item, branchless, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<ItemType>& keys, const size_t start, const size_t end, Table* table) {
    table->AddAll(keys, start, end);
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return (0 == table->Contain(key));
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = SuccinctCountingBlockedBloomFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Contain(key);
  }
};

template <typename ItemType, size_t bits_per_item, typename HashFamily>
struct FilterAPI<SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily>> {
  using Table = SuccinctCountingBlockedBloomRankFilter<ItemType, bits_per_item, HashFamily>;
  static Table ConstructFromAddCount(size_t add_count) { return Table(add_count); }
  static void Add(uint64_t key, Table* table) {
    table->Add(key);
  }
  static void AddAll(const vector<uint64_t>& keys, const size_t start, const size_t end, Table* table) {
    for(size_t i = start; i < end; i++) { Add(keys[i],table); }
  }
  static void Remove(uint64_t key, Table * table) {
    table->Remove(key);
  }
  CONTAIN_ATTRIBUTES static bool Contain(uint64_t key, const Table * table) {
    return table->Contain(key);
  }
};


#endif