# xorfilter_cpp
Xor Filters: Faster and Smaller Than Bloom and Cuckoo Filters (C++)

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

See src/xorfilter.h. This single header depends on src/hashutil.h.

## Credit

The code is derived from https://github.com/efficient/cuckoofilter by Bin Fan et al.