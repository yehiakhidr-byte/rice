#!/bin/bash
cd /home/yehia/sensitivity
for f in sens_D*.bim; do
  echo "$f: $(wc -l < $f)"
done
echo "---"
# Also check the raw QC file
find /home/yehia -name "rice_QC.bim" 2>/dev/null
for f in rice_pruned?.log; do
  grep "variants loaded" $f 2>/dev/null
done