#!/bin/bash
# HDRA Post-processing: FST, diversity, LD decay
# Run after hdra_cont.sh completes
set -e
cd /home/yehia/hdra_rdp1

echo "=== HDRA Post-Pipeline Analysis ==="
echo "Scenario,SNPs,Mean_FST" > hdra_fst_summary.csv

# Read subpop map
awk '{print $1, $1, $2}' subpop_map.txt > hdra_pop.cluster

for tag in hdra_S_none hdra_S_010 hdra_S_005 \
           hdra_M_none hdra_M_001 hdra_M_003 hdra_M_005 hdra_M_010 \
           hdra_G_none hdra_G_020 hdra_G_010 hdra_G_005 hdra_G_002 \
           hdra_L_none hdra_L_08 hdra_L_05 hdra_L_02 hdra_L_100; do
  
  # Determine input file
  if [ -f "${tag}.pruned.bed" ]; then
    bfile="${tag}.pruned"
  elif [ -f "${tag}.bed" ]; then
    bfile="${tag}"
  else
    echo "  $tag: no data, skipping"
    continue
  fi
  
  n_snps=$(wc -l < ${bfile}.bim)
  n_samps=$(wc -l < ${bfile}.fam)
  echo "  $tag: $n_samps x $n_snps SNPs"
  
  # FST
  if [ ! -f "${tag}_fst.fst" ]; then
    echo "    Computing FST..."
    plink1.9 --bfile $bfile --fst --within hdra_pop.cluster --out ${tag}_fst 2>&1 | tail -2
    if [ -f "${tag}_fst.fst" ]; then
      fst_val=$(awk 'NR>1 && !/NaN/ {sum+=$NF; n++} END {if(n>0) printf "%.4f", sum/n; else print "NA"}' ${tag}_fst.fst)
      echo "$tag,$n_snps,$fst_val" >> hdra_fst_summary.csv
    fi
  else
    echo "    FST already done"
  fi
  
  # LD decay for L scenarios (and baseline)
  if [[ "$tag" == hdra_L_* ]] || [[ "$tag" == "hdra_S_010" ]]; then
    if [ ! -f "${tag}_ld.ld" ]; then
      echo "    Computing LD decay..."
      n_snps2=$(wc -l < ${bfile}.bim)
      if [ "$n_snps2" -gt 50000 ]; then
        # Thin to every 10th SNP
        awk 'NR % 10 == 0 {print $2}' ${bfile}.bim > ${tag}_thin.snps
        plink1.9 --bfile $bfile --extract ${tag}_thin.snps --make-bed --out ${tag}_thin 2>&1 | tail -1
        ld_bfile="${tag}_thin"
      else
        ld_bfile="${bfile}"
      fi
      # LD on chromosome 1 only for speed
      echo "    LD on chr1 ($(awk '$1==1' ${ld_bfile}.bim | wc -l) SNPs)..."
      plink1.9 --bfile $ld_bfile --chr 1 --r2 --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0.0 --out ${tag}_ld 2>&1 | tail -3
    else
      echo "    LD already computed"
    fi
  fi
done

echo ""
echo "=== FST Summary ==="
cat hdra_fst_summary.csv
echo ""
echo "=== Post-processing Complete ==="