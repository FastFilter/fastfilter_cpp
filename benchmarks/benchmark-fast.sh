#!/bin/sh
# run the benchmark multiple times with all important algorithms
# for algorithm ids and other parameters, see
# bulk-insert-and-query.cc
#
# rnd: random number generators to use
for rnd in `seq -1 -1`; do
  # m: number of entries (in millions)
  for m in `seq 100 100`; do
    # test: test id
    for test in `seq 1 10`; do
      sleep 5;
      # default algorithms
      ./bulk-insert-and-query.exe ${m}000000;
    done;
  done;
done > benchmark-fast-results.txt 2>&1
