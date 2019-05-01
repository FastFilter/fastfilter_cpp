#!/bin/sh
for rnd in `seq -1 -1`; do
  for alg in 0 1 2 3 4 10 11 12 13 20 40 41 42 43 44 45 46 47 48 50 51 52 60 61 62 100; do
    for m in `seq 10 90 100`; do
      for test in `seq 1 3`; do
        now=$(date +"%T");
        echo ${now} alg ${alg} size ${m} ${rnd};
        ./bulk-insert-and-query.exe ${m}000000 ${alg} ${rnd};
      done;
    done;
  done;
done > benchmark-results.txt
