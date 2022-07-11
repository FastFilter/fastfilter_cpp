#ifndef HASHUTIL_H_
#define HASHUTIL_H_

#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>

#include <string>
#include <fstream>

#include <random>

namespace hashing {

inline size_t sysrandom(void *dst, size_t dstlen) {
        char *buffer = reinterpret_cast<char *>(dst);
        std::ifstream stream("/dev/urandom", std::ios_base::binary | std::ios_base::in);
        stream.read(buffer, dstlen);

        return dstlen;
    }

    // See Martin Dietzfelbinger, "Universal hashing and k-wise independent random
    // variables via integer arithmetic without primes".
    class TwoIndependentMultiplyShift {
        unsigned __int128 multiply_, add_;

    public:
        TwoIndependentMultiplyShift() {
            std::uint_least64_t seed;
            sysrandom(&seed, sizeof(seed));
            std::mt19937_64 rng(seed);
            std::uniform_int_distribution<uint64_t> dist(0, UINT64_MAX);

            for (auto v : {&multiply_, &add_}) {
                for (size_t i = 0; i < 8; i++) {
                    unsigned __int128 hi = dist(rng);
                    unsigned __int128 lo = dist(rng);
                    *v ^= (hi << 64) | lo;
                }
            }
        }
        
        inline uint64_t operator()(uint64_t key) const {
            return (add_ + multiply_ * static_cast<decltype(multiply_)>(key)) >> 64;
        }
    };


class SimpleMixSplit {

 public:
  uint64_t seed;
  SimpleMixSplit() {
    ::std::random_device random;
    seed = random();
    seed <<= 32;
    seed |= random();
  }

  inline static uint64_t murmur64(uint64_t h) {
    h ^= h >> 33;
    h *= UINT64_C(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h *= UINT64_C(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    return h;
  }

  inline uint64_t operator()(uint64_t key) const {
    return murmur64(key + seed);
  }
};

}

#endif  // CUCKOO_FILTER_HASHUTIL_H_
