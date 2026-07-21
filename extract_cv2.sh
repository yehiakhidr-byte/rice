#!/bin/bash
cd /root
for K in 1 2 3 4 5 6 7 8 9 10; do
  echo -n "K=$K: "
  grep "CV error" log2_${K}.out
done