// Generating random data

#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <random>
#include <stdexcept>
#include <vector>

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
  return result;
}



/**********************************************
** WARNING: MixIn and MixInFast can generate duplicates!!!! Use SimpleMixIn and SlowishSimpleMixIn
** instead.
*******************************************/
// Using two pointer ranges for sequences x and y, create a vector clone of x but for
// y_probability y's mixed in.
template <typename T>
::std::vector<T> MixIn(const T* x_begin, const T* x_end, const T* y_begin, const T* y_end,
    double y_probability) {
  const size_t x_size = x_end - x_begin, y_size = y_end - y_begin;
  if (y_size > (1ull << 32)) throw ::std::length_error("y is too long");
  ::std::vector<T> result(x_begin, x_end);
  ::std::random_device random;
  auto genrand = [&random, y_size]() {
    return (static_cast<size_t>(random()) * y_size) >> 32;
  };
  for (size_t i = 0; i < y_probability * x_size; ++i) {
    result[i] = *(y_begin + genrand());
  }
  ::std::shuffle(result.begin(), result.end(), random);
  return result;
}

template <typename T>
::std::vector<T> MixInFast(const T* x_begin, const T* x_end, const T* y_begin, const T* y_end,
    double y_probability, uint64_t start) {
  const size_t x_size = x_end - x_begin, y_size = y_end - y_begin;
  if (y_size > (1ull << 32)) throw ::std::length_error("y is too long");
  ::std::vector<T> result(x_begin, x_end);
  uint64_t index = start;
  auto genrand = [&index, y_size]() {
    // mix64
    uint64_t x = index++;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9L;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebL;
    x = x ^ (x >> 31);
    return static_cast<size_t>((((uint32_t) x) * (uint64_t) y_size) >> 32);
  };
  for (size_t i = 0; i < y_probability * x_size; ++i) {
    result[i] = *(y_begin + genrand());
  }
  ::std::shuffle(result.begin(), result.end(), std::default_random_engine(start));
  return result;
}
