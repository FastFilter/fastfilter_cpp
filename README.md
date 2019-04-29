# fastfilter_cpp

Fast Filter: Fast approximate membership filter implementations (C++)

This is a research library currently. Some of the source code is not fully tested, and not well documented. 

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
                                                    find    find    find    find    find                        optimal  wasted million
                                     add  remove      0%     25%     50%     75%    100%        ε  bits/item  bits/item   space    keys

add    cycles: 325.5/key, instructions: (303.2/key, 0.93/cycle) cache misses: 12.41/key branch misses: 1.17/key
0.00%  cycles:  81.7/key, instructions: ( 48.0/key, 0.59/cycle) cache misses:  3.06/key branch misses: 0.00/key
0.25%  cycles:  81.8/key, instructions: ( 48.0/key, 0.59/cycle) cache misses:  3.06/key branch misses: 0.00/key
0.50%  cycles:  81.8/key, instructions: ( 48.0/key, 0.59/cycle) cache misses:  3.06/key branch misses: 0.00/key
0.75%  cycles:  82.0/key, instructions: ( 48.0/key, 0.59/cycle) cache misses:  3.06/key branch misses: 0.00/key
1.00%  cycles:  81.9/key, instructions: ( 48.0/key, 0.59/cycle) cache misses:  3.06/key branch misses: 0.00/key
                            Xor8  106.79    0.00   25.92   25.88   25.86   25.94   25.98  0.3892%       9.84       8.01   22.9%    10.0
                            
... # many more lines omitted
```

The `add` lines preceding the name of each algorithm gives you information regarding the construction time whereas
the other five lines give you information regarding the queries where a given percentage of elements are present
in the set. We use Linux performance counters to measure instructions, cache misses and branch misses.

As part of the benchmark, we check the correctness of the implementation.


## Where is your code?

The filter implementations are in `src`, most are single header files and depend on `src/hashutil.h`:

* src/bloom.h
* src/xorfilter.h

## Credit

The cuckoo filter and the benchmark are derived from https://github.com/efficient/cuckoofilter by Bin Fan et al.
