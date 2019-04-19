# fastfilter_cpp
Fast Filter: Fast approximate membership filter implementations (C++)

This is a research library, developers might want to consider our [Header-only Xor Filter library in C](https://github.com/FastFilter/xor_singleheader/).

## Prerequisites

- A  C++11 compiler such as GNU G++ or LLVM Clang++
- Make 



## Usage

```
cd benchmarks
make
./bulk-insert-and-query.exe 10000000
```


## Where is your code?

The filter implementations are in `src`, most are single header files and depend on `src/hashutil.h`:

* src/bloom.h
* src/xorfilter.h

## Credit

The cuckoo filter and the benchmark are derived from https://github.com/efficient/cuckoofilter by Bin Fan et al.
