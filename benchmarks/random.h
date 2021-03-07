// Generating random data

#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <random>
#include <stdexcept>
#include <vector>

// this can be atrociously slow
template <class RNG = ::std::random_device>
::std::vector<::std::uint64_t> GenerateRandom64(::std::size_t count) {
  ::std::vector<::std::uint64_t> result(count);
  RNG random;
  // To generate random keys to lookup, this uses ::std::random_device which is slower but
  // stronger than some other pseudo-random alternatives. The reason is that some of these
  // alternatives (like libstdc++'s ::std::default_random, which is a linear congruential
  // generator) behave non-randomly under some hash families like Dietzfelbinger's
  // multiply-shift.
  auto genrand = [&random]() {
    return random() + (static_cast<::std::uint64_t>(random()) << 32);
  };
  ::std::generate(result.begin(), result.end(), ::std::ref(genrand));
  return result;
}

::std::vector<::std::uint64_t> GenerateRandom64Fast(::std::size_t count, uint64_t start) {
  ::std::vector<::std::uint64_t> result(count);
  uint64_t index = start;
  auto genrand = [&index]() {
    // mix64
    uint64_t x = index++;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9L;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebL;
    x = x ^ (x >> 31);
    return x;
  };
  ::std::generate(result.begin(), result.end(), ::std::ref(genrand));
  return result;
}



static inline uint64_t biased_random_bounded(uint64_t range, __uint128_t * seed) {
  __uint128_t random64bit, multiresult;
  *seed *= UINT64_C(0xda942042e4dd58b5);
  random64bit = *seed >> 64;
  multiresult = random64bit * range;
  return multiresult >> 64; // [0, range)
}

static inline uint64_t random_bounded(uint64_t range, __uint128_t * seed) {
  __uint128_t random64bit, multiresult;
  uint64_t leftover;
  uint64_t threshold;
  *seed *= UINT64_C(0xda942042e4dd58b5);
  random64bit = *seed >> 64;
  multiresult = random64bit * range;
  leftover = (uint64_t)multiresult;
  if (leftover < range) {
    threshold = -range % range;
    while (leftover < threshold) {
      *seed *= UINT64_C(0xda942042e4dd58b5);
      random64bit = *seed >> 64;
      multiresult = random64bit * range;
      leftover = (uint64_t)multiresult;
    }
  }
  return multiresult >> 64; // [0, range)
}

template <typename T>
// pick capacity elements form x_begin, x_end, write them to storage
void reservoirsampling(T *storage, uint32_t capacity, const T* x_begin, const T* x_end, __uint128_t * seed) {
  if(capacity == 0) return;
  size_t size = x_end - x_begin;
  if(size < capacity) {
    throw ::std::logic_error("I cannot sample the requested number. This is not going to end well.");
  }
  size_t i;
  for (i = 0; i < capacity; i++) {
    storage[i] = x_begin[i];
  }
  while(i < size) {
    size_t nextpos =
        biased_random_bounded(i, seed);
    if(nextpos < capacity) {
      storage[nextpos] = x_begin[i];
    }
    i++;
  }
}


// good old Fisher-Yates shuffle, shuffling an array of integers, without
// division
template <typename T>
void fast_shuffle(T *storage, uint64_t size, __uint128_t* seed) {
  uint64_t i;
  for (i = size; i > 1; i--) {
    uint64_t nextpos = random_bounded(i, seed);
    T tmp = storage[i - 1];   // likely in cache
    T val = storage[nextpos]; // could be costly
    storage[i - 1] = val;
    storage[nextpos] = tmp; // you might have to read this store later
  }
}

// Using two pointer ranges for sequences x and y, create a vector clone of x but for
// y_probability y's mixed in.
template <typename T>
::std::vector<T> DuplicateFreeMixIn(const T* x_begin, const T* x_end, const T* y_begin, const T* y_end,
    double y_probability, uint64_t start) {
  const size_t x_size = x_end - x_begin; //, y_size = y_end - y_begin;
  ::std::vector<T> result;
  result.resize(x_size);
  __uint128_t seed = start;
  size_t howmanyy = round(x_size * y_probability);
  size_t howmanyx = x_size - howmanyy;
  reservoirsampling(result.data(), howmanyx,  x_begin, x_end, &seed);
  reservoirsampling(result.data() + howmanyx, howmanyy,  y_begin, y_end, &seed);
  if((y_probability != 0.0) && (y_probability != 1.0)) { fast_shuffle(result.data(), result.size(), &seed); }
  return result;
}


