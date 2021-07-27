
#include "filterapi.h"
#include "random.h"
#include "timing.h"

struct simple_stat {
  double avg_bits_per_item;
  double avg_time_per_item;
};

template <typename Table = xorbinaryfusefilter_lowmem::XorBinaryFuseFilter<
              uint64_t, uint8_t>>
simple_stat
one_measure_benchmark_construction(const std::vector<uint64_t> &source,
                                   size_t add_count) {
  if (add_count >= source.size()) {
    throw std::runtime_error("insufficient capacity");
  }
  Table filter = FilterAPI<Table>::ConstructFromAddCount(add_count);
  size_t total_time{0};
  double total_bits_per_item{0};
  size_t iterations{0};
  // We run each test so that they take at least 0.1 second.
  while (total_time < 100'000'000) {
    auto start_time = NowNanos();
    FilterAPI<Table>::AddAll(source, 0, add_count, &filter);
    auto time = NowNanos() - start_time;
    double bits_per_item = double(filter.SizeInBytes()) * 8.0 / add_count;
    total_time += time;
    total_bits_per_item += bits_per_item;
    iterations++;
  }
  double avg_bits_per_item = total_bits_per_item / iterations;
  double avg_time_per_item = total_time / (iterations * add_count);
  return {avg_bits_per_item, avg_time_per_item};
}

void stream_benchmark(size_t start_size, size_t max_size,
                      size_t iteration = 100) {
  std::vector<uint64_t> source = GenerateRandom64Fast(max_size, rand());
  std::cout << "n, binbits, bintime, bin4bits, bin4time, xorbits, xortime, bloombits, bloomtime" << std::endl;

  double gap = exp(log(max_size / start_size) / iteration);
  for (double n = start_size; n <= max_size; n *= gap) {
    std::cout << size_t(n) << ",\t";
    auto s = one_measure_benchmark_construction<>(source, size_t(n));
    std::cout << s.avg_bits_per_item << "," << s.avg_time_per_item;
    std::cout.flush();
    std::cout << ",\t";
    s = one_measure_benchmark_construction<xorbinaryfusefilter_lowmem4wise::XorBinaryFuseFilter<uint64_t, uint8_t>>(source, size_t(n));
    std::cout << s.avg_bits_per_item << "," << s.avg_time_per_item;
    std::cout.flush();
    std::cout << ",\t";
    s = one_measure_benchmark_construction<
        XorFilter<uint64_t, uint8_t, SimpleMixSplit>>(source, size_t(n));
    std::cout << s.avg_bits_per_item << "," << s.avg_time_per_item;
    std::cout << ",\t";
    s = one_measure_benchmark_construction<
        BloomFilter<uint64_t, 12, false, SimpleMixSplit>>(source, size_t(n));
    std::cout << s.avg_bits_per_item << "," << s.avg_time_per_item;
    std::cout << std::endl;
  }
}
int main() {
  stream_benchmark(100'000, 100'000'000);
  return EXIT_SUCCESS;
}