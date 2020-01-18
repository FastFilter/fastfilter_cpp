#!/bin/sh
# runs the benchmark program once with all important algorithms,
# in order to get the space overhead for certain sizes
# (not to get timing data)
# for algorithm ids and other parameters, see
# bulk-insert-and-query.cc
#
# rnd: random number generators to use
for rnd in `seq -1 -1`; do
  # alg: algorithms to test (not used currently)
  for alg in 0; do
    # m: number of entries
    for m in `seq 10 1 100`; do
      # test: test id
      for test in `seq 1 1`; do
        now=$(date +"%T");
        echo ${now} alg ${alg} size ${m} ${rnd};
        ./bulk-insert-and-query.exe ${m}00000 "0,2,3,4,11,12,13,15,16,17,20,30,40,41,42,51,80" ${rnd};
        #./bulk-insert-and-query.exe ${m}00000 "13,15,16,17" ${rnd};
      done;
    done;
  done;
done> spaceUsage-results.txt 2>&1
