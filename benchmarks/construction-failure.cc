
#include "filterapi.h"
#include "random.h"
#include "timing.h"


template <typename Table = xorbinaryfusefilter_lowmem::XorBinaryFuseFilter<
              uint64_t, uint8_t>>
std::pair<double, double>
one_measure_benchmark_construction(const std::vector<uint64_t> &source,
                                   size_t add_count) {
  if (add_count >= source.size()) {
    throw std::runtime_error("insufficient capacity");
  }
  Table filter = FilterAPI<Table>::ConstructFromAddCount(add_count);
  size_t total{1'000};
  if(add_count > 50'000) { total = 100; }
  if(add_count > 500'000) { total = 10; }
  size_t iterations{0};
  size_t failures{0};
  double bits_per_item{0.0};
  // We run each test so that they take at least 0.1 second.
  while (iterations < total) {
    FilterAPI<Table>::AddAll(source, 0, add_count, &filter);
    failures += filter.hashIndex;
    bits_per_item += double(filter.SizeInBytes()) * 8.0 / add_count;

    iterations++;
  }
  size_t total_attempts = total + failures;
  return make_pair(double(failures) / double(total_attempts), bits_per_item / total);
}

void stream_construction_benchmark(size_t start_size, size_t max_size,
                      size_t iteration = 100) {
  std::vector<uint64_t> source = GenerateRandom64Fast(max_size, rand());
  std::cout << "n, binfail, binbits, bin4fail, bin4bits,  xorfail, xorbits" << std::endl;
  std::cout << std::setprecision(5) << std::fixed;

  double gap = exp(log(max_size / start_size) / iteration);
  for (double n = start_size; n <= max_size; n *= gap) {
    std::cout << size_t(n) << ",\t";
    std::pair<double,double> failure_proba = one_measure_benchmark_construction<>(source, size_t(n));
    std::cout << failure_proba.first << "," << failure_proba.second;
    std::cout.flush();
    std::cout << ",\t";
    failure_proba = one_measure_benchmark_construction<xorbinaryfusefilter_lowmem4wise::XorBinaryFuseFilter<uint64_t, uint8_t>>(source, size_t(n));
    std::cout << failure_proba.first << "," << failure_proba.second;
    std::cout.flush();
    std::cout << ",\t";
    failure_proba = one_measure_benchmark_construction<XorFilter<uint64_t, uint8_t, SimpleMixSplit>>(source, size_t(n));
    std::cout << failure_proba.first << "," << failure_proba.second;
    std::cout << std::endl;
  }
}
int main() {
  stream_construction_benchmark(100, 100'000, 200);
  return EXIT_SUCCESS;
}