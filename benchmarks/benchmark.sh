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
    for test in `seq 1 20`; do
      # alg: algorithms to test
      for alg in 0 2 11 12 44 45 51 116 117 118 119 1086 1156 3086 3156; do
        now=$(date +"%T");
        echo ${test} ${now} alg ${alg} size ${m} ${rnd};
        sleep 10;
        ./bulk-insert-and-query.exe ${m}000000 ${alg} ${rnd};
      done;
    done;
  done;
done > benchmark-results.txt 2>&1
