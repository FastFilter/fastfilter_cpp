#ifndef HASHUTIL_H_
#define HASHUTIL_H_

#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>

#include <string>

#include <random>

namespace hashing {
// See Martin Dietzfelbinger, "Universal hashing and k-wise independent random
// variables via integer arithmetic without primes".

class TwoIndependentMultiplyShift {
  unsigned __int128 multiply_, add_;
 public:
  TwoIndependentMultiplyShift() {
    ::std::random_device random;
    for (auto v : {&multiply_, &add_}) {
      *v = random();
      for (int i = 1; i <= 4; ++i) {
        *v = *v << 32;
        *v |= random();
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