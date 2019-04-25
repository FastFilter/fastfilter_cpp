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
#ifndef _FIXED_POINT_H
#define _FIXED_POINT_H

namespace CompressedCuckoo{
  // This class is used for getting around the restriction that template 
  // parameters need to be integer types.  I want a compile-time floating 
  // point template parameter, so I'm using this.
  typedef __uint128_t SerializedFixedPoint;
  struct FixedPoint{
    uint64_t _numerator;
    uint64_t _denominator;
    constexpr FixedPoint(uint64_t numerator, uint64_t denominator) : 
      _numerator(numerator), _denominator(denominator){}
    // Specialized constructor for values between 0 and 1
    explicit constexpr FixedPoint(double fp_representation) :
      _numerator(fp_representation * 0x8000000000000000llu), 
      _denominator(0x8000000000000000llu){}  
    explicit constexpr FixedPoint(SerializedFixedPoint serialized_fixed_point) :
      _numerator(serialized_fixed_point & (0xffffffffffffffffllu)), 
      _denominator(serialized_fixed_point >> 64){}
    constexpr double to_double() const{
      return static_cast<double>(_numerator) / _denominator;
    }
    constexpr float to_float() const{
      return static_cast<float>(_numerator) / _denominator;
    }
    constexpr SerializedFixedPoint serialize() const{
      return _numerator + (static_cast<SerializedFixedPoint>(_denominator) 
        << 64);
    }
  };
}

#endif
