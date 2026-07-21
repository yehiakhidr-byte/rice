#!/bin/bash
cd /home/yehia/sensitivity
for f in sens_D1_S2 sens_D1_S3 sens_D2_M0 sens_D2_M1 sens_D2_M2 sens_D2_M4 sens_D3_G002 sens_D3_G010 sens_D3_G020 sens_D3_Gnone sens_D4_LD05 sens_D4_LD08 sens_D4_LD100 sens_D4_LDnone rice_pruned rice_pruned2 rice_pruned3; do
  if [ -f "${f}.bed" ] && [ ! -f "${f}_dist.txt" ]; then
    echo "Dist: $f"
    plink1.9 --bfile "$f" --distance square 1-ibs --out "${f}_dist" 2>&1 | tail -1
  fi
done
echo "Done"
ls -la /home/yehia/sensitivity/*_dist.txt 2>/dev/null
