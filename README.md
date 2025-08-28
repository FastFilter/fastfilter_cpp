# fastfilter_cpp
[![Ubuntu 22.04 CI (GCC 11)](https://github.com/FastFilter/fastfilter_cpp/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/FastFilter/fastfilter_cpp/actions/workflows/ubuntu.yml)

Fast Filter: Fast approximate membership filter implementations (C++)

This is a research library currently. It is not meant for production use.

Developers might want to consider our [Header-only Xor Filter library in C](https://github.com/FastFilter/xor_singleheader/).


## Prerequisites

- A  C++11 compiler such as GNU G++ or LLVM Clang++
- Make

Expectations:

- Though it should be possible to run this benchmark on any operating system, we expect Linux and use its performance counters to measure performance.
- We expect an x64 processor with AVX2 support though most filters work on any processor, if you compile on a machine that does not support AVX2 instructions, the corresponding filters that depend on AVX2 will be disabled.

## Usage

Make sure to select the right GNU GCC compiler (e.g., via `export export CXX=g++-8`).
You may want to disable hyperthreading and adjust page sizes. Run the benchmark
on a quiet machine.


```
git clone https://github.com/FastFilter/fastfilter_cpp.git
cd fastfilter_cpp
cd benchmarks
make
# there may be compiler warnings at this point, we compile with '-Wall'
./bulk-insert-and-query.exe 10000000
# collect the output (it is quite verbose)
./bulk-insert-and-query.exe 100000000
```

Your results will depend on the hardware, on the compiler and how the system is configured. A sample output is as follows:

```
$ ./bulk-insert-and-query.exe 10000000
./bulk-insert-and-query.exe 10000000
                                                    find    find    find    find    find  1*add+                       optimal  wasted million
                                     add  remove      0%     25%     50%     75%    100%  3*find      ε%  bits/item  bits/item  space%    keys

add    cycles: 351.3/key, instructions: (332.4/key, 0.95/cycle) cache misses: 13.99/key branch misses: 1.23/key
0.00%  cycles:  87.8/key, instructions: ( 48.0/key, 0.55/cycle) cache misses:  2.89/key branch misses: 0.00/key
0.25%  cycles:  87.8/key, instructions: ( 48.0/key, 0.55/cycle) cache misses:  2.89/key branch misses: 0.00/key
0.50%  cycles:  87.8/key, instructions: ( 48.0/key, 0.55/cycle) cache misses:  2.90/key branch misses: 0.00/key
0.75%  cycles:  87.9/key, instructions: ( 48.0/key, 0.55/cycle) cache misses:  2.90/key branch misses: 0.00/key
1.00%  cycles:  87.8/key, instructions: ( 48.0/key, 0.55/cycle) cache misses:  2.89/key branch misses: 0.00/key
                            Xor8  106.59    0.00   23.85   23.85   23.85   23.88   23.83  178.15  0.3908       9.84       8.00    23.0  10.000
... # many more lines omitted
```

The `add` lines preceding the name of each algorithm gives you information regarding the construction time whereas
the other five lines give you information regarding the queries where a given percentage of elements are present
in the set. We use Linux performance counters to measure instructions, cache misses and branch misses.

As part of the benchmark, we check the correctness of the implementation.

## Benchmarking

A simple way to run the benchmark for the most important algorithms
is to run `./bulk-insert-and-query.exe <number of entries>`.
To get a low error, it is best run on a Linux machine that is not otherwise in use.
You might want to run the tests multiple times to verify it is working properly.
The steps to run the tests are:

    git clone https://github.com/FastFilter/fastfilter_cpp.git
    cd fastfilter_cpp/benchmarks
    make clean ; make
    # 100 million entries
    ./bulk-insert-and-query.exe 100000000

## Where is your code?

The filter implementations are in `src/<type>/`. Most implementations depend on `src/hashutil.h`. Examples:

* src/bloom/bloom.h
* src/xorfilter/xorfilter.h
* src/xorfilter/3wise_xor_binary_fuse_filter_lowmem.h


## References

- [Binary Fuse Filters: Fast and Smaller Than Xor Filters](http://arxiv.org/abs/2201.01174), Journal of Experimental Algorithmics 27, 2022.
- [Xor Filters: Faster and Smaller Than Bloom and Cuckoo Filters](https://arxiv.org/abs/1912.08258), Journal of Experimental Algorithmics 25 (1), 2020
- [Prefix Filter: Practically and Theoretically Better Than Bloom](https://arxiv.org/abs/2203.17139), PVLDB 15(7), 2022.

## Credit

The cuckoo filter and the benchmark are derived from https://github.com/efficient/cuckoofilter by Bin Fan et al.
The SIMD blocked Bloom filter is from https://github.com/apache/impala (via the cuckoo filter).
The Morton filter is from https://github.com/AMDComputeLibraries/morton_filter.
The counting quotient filter (CQF) is from https://github.com/splatlab/cqf.
The vector quotient filter is from https://github.com/splatlab/vqf. 
The ribbon filters are from https://github.com/pdillinger/fastfilter_cpp.
The prefix filter is from https://github.com/TomerEven/Prefix-Filter.


# Implementations of xor and binary fuse filters in other programming languages

* [C](https://github.com/FastFilter/xor_singleheader)
* [C99](https://github.com/skeeto/xf8)
* [Erlang](https://github.com/mpope9/exor_filter)
* [Go](https://github.com/FastFilter/xorfilter)
* [Java](https://github.com/FastFilter/fastfilter_java)
* [Python](https://github.com/GreyDireWolf/pyxorfilter)
* Rust: [1](https://github.com/bnclabs/xorfilter), [2](https://github.com/codri/xorfilter-rs), [3](https://github.com/Polochon-street/rustxorfilter)
* [Zig](https://github.com/hexops/xorfilter)
* [Julia](https://github.com/JokingHero/FastFilter.jl)
* [C#](https://github.com/jonmat/FastIndex)
