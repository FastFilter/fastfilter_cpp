#include "bloom.h"
#include <random>
#include <iostream>

bool check_bloom() {
    std::cout << "testing Bloom" << std::endl;
	bloomfilter::BloomFilter<uint64_t, 18, false> bf(65536);
	for (uint64_t i = 0; i < 65536; i++) {
		bf.Add(i);
	}
    std::mt19937 gen;
    std::discrete_distribution<> bytes_count;
    std::uniform_int_distribution<uint64_t> val_64bit{0x0, 0xffffffffffffffff};
	uint64_t cnt = 0;
    uint64_t count = 50000000;
	for (uint64_t i = 0; i < count; i++)
	{
		int ans = bf.Contain(val_64bit(gen));
		if (ans != 1) { cnt++; }
	}
    double rate = double(cnt)/count;
	std::cout << double(cnt)/count << std::endl;
    return rate < 0.00018;
}

int main() {
    return check_bloom() ? EXIT_SUCCESS : EXIT_FAILURE;
}