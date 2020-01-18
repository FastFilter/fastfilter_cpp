#!/bin/sh
# run the benchmark multiple times with all important algorithms
# for algorithm ids and other parameters, see
# bulk-insert-and-query.cc
#
# rnd: random number generators to use
for rnd in `seq -1 -1`; do
  # alg: algorithms to test
  for alg in 0 2 3 4 11 12 13 15 16 17 20 30 40 41 42 51 80 100; do
    # m: number of entries
    for m in `seq 10 90 100`; do
      # test: test id
      for test in `seq 1 3`; do
        now=$(date +"%T");
        echo ${now} alg ${alg} size ${m} ${rnd};
        ./bulk-insert-and-query.exe ${m}000000 ${alg} ${rnd};
      done;
    done;
  done;
done > benchmark-results.txt 2>&1
