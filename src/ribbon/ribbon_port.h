//  Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).

#pragma once

#define CACHE_LINE_SIZE 64U

namespace port {

// FIXME
constexpr bool kLittleEndian = true;

#define PREFETCH(addr, rw, locality) __builtin_prefetch(addr, rw, locality)

#define HAVE_AVX2 __AVX2__

}  // namespace port
