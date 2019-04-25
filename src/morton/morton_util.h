/*
Copyright (c) 2019 Advanced Micro Devices, Inc.
 
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Author: Alex D. Breslow 
        Advanced Micro Devices, Inc.
        AMD Research
*/
#ifndef _MORTON_UTIL_H
#define _MORTON_UTIL_H

namespace CompressedCuckoo{
  // Computes the ceiling on input/divisor and then adds 1 if necessary to make 
  // the result even
  template<class T>
  T divide_round_up_even(T input, T divisor){
    T ret = (input + divisor - 1) / divisor;
    return ret & 1 ? ret + 1 : ret;
  }

  // total_slots should be the initial estimate for the number of total logical
  // slots in the filter, not physical slots
  template<class T>
  T determine_total_buckets(T slots_per_bucket, T total_slots, 
    T buckets_per_block){
    // Alternate implementation that's slightly different in that it doesn't 
    // round up to the next whole block which would produce an even number of 
    // buckets.  Instead, it just rounds up by 1 if the number of buckets would 
    // be odd.
    // Provisional total slots
    //total_slots = ((total_slots + slots_per_bucket - 1) / slots_per_bucket) * 
    //  slots_per_bucket;
    //return divide_round_up_even(total_slots, slots_per_bucket);

    T logical_slots_per_block = slots_per_bucket * buckets_per_block;
    T total_blocks = (total_slots + logical_slots_per_block - 1) / 
      (logical_slots_per_block);
    T total_buckets = total_blocks * buckets_per_block;
    total_buckets = (total_buckets & 1) ? total_buckets + buckets_per_block : 
      total_buckets; 
    return total_buckets;
  }
}

#endif // End of file guards
