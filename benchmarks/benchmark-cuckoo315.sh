#!/bin/sh
# run the benchmark multiple times with all important algorithms
# for algorithm ids and other parameters, see
# bulk-insert-and-query.cc
#
# rnd: random number generators to use
for rnd in `seq -1 -1`; do
  # alg: algorithms to test
  for alg in 11 12 13 15 16 17; do
    # m: number of entries
    for m in 315; do
      # test: test id
      for test in `seq 1 3`; do
        now=$(date +"%T");
        echo ${now} alg ${alg} size ${m} ${rnd};
        ./bulk-insert-and-query.exe ${m}00000 ${alg} ${rnd};
      done;
    done;
  done;
done > benchmark-results.txt 2>&1
