#!/bin/bash
for K in $(seq 1 12); do
  best=$(cat /root/admixture_results/K${K}/best.txt)
  ll=$(echo $best | awk '{print $1}')
  cv=$(echo $best | awk '{print $2}')
  printf "K=%-3d CV=%.5f  LL=%.2f\n" $K $cv $ll
done
