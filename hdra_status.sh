#!/bin/bash
cd /home/yehia/hdra_rdp1
echo "=== Processes ==="
ps aux | grep -E 'hdra_cont|admixture' | grep -v grep
echo "=== ADMIXTURE Progress ==="
for s in hdra_S_none hdra_S_010 hdra_S_005 hdra_M_none hdra_M_001 hdra_M_003 hdra_M_005 hdra_M_010 hdra_G_none hdra_G_020 hdra_G_010 hdra_G_005 hdra_G_002 hdra_L_none hdra_L_08 hdra_L_05 hdra_L_02 hdra_L_100; do
  nq=$(ls ${s}_K*.Q 2>/dev/null | wc -l)
  [ "$nq" -gt 0 ] && echo "  $s: ${nq}/10"
done
echo "=== Status ==="
cat p.csv 2>/dev/null
echo "=== Stats ==="
cat s.csv 2>/dev/null