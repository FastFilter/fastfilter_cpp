#ifndef XOR_BINARY_FUSE_FILTER_XOR_FILTER_H_
#define XOR_BINARY_FUSE_FILTER_XOR_FILTER_H_

#include <algorithm>
#include <assert.h>

#include "hashutil.h"
#include "timing.h"

using namespace std;
using namespace hashing;

size_t calculateSegmentLength(size_t arity, size_t size) {
  size_t segmentLength;
  if (arity == 3) {
    // We deliberately divide a log by a log so that the reader does not have
    // to ask about the basis of the log.
    segmentLength = 1L << (int)floor(log(size) / log(3.33) + 2.25);
  } else if (arity == 4) {
    segmentLength = 1L << (int)floor(log(size) / log(2.91) - 0.5);
  } else {
    // not supported
    segmentLength = 65536;
  }
  return segmentLength;
}

double calculateSizeFactor(size_t arity, size_t size) {
  double sizeFactor;
  if (arity == 3) {
    sizeFactor = fmax(1.125, 0.875 + 0.25 * log(1000000) / log(size));
  } else if (arity == 4) {
    sizeFactor = fmax(1.075, 0.77 + 0.305 * log(600000) / log(size));
  } else {
    // not supported
    sizeFactor = 2.0;
  }
  return sizeFactor;
}

#include "3wise_xor_binary_fuse_filter_naive.h"
#include "3wise_xor_binary_fuse_filter_lowmem.h"

#include "4wise_xor_binary_fuse_filter_naive.h"
#include "4wise_xor_binary_fuse_filter_lowmem.h"

#endif // XOR_BINARY_FUSE_FILTER_XOR_FILTER_H_
