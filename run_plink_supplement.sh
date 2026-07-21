#!/bin/bash
# Run PLINK supplementary analyses for sensitivity datasets
# PCA, freq, HWE for all datasets; LD decay for LD dimension

SENS_DIR="/home/yehia/sensitivity"
DATASETS="sens_D1_S2 sens_D1_S3 sens_D2_M0 sens_D2_M1 sens_D2_M2 sens_D2_M4 sens_D3_G002 sens_D3_G010 sens_D3_G020 sens_D3_Gnone sens_D4_LD05 sens_D4_LD08 sens_D4_LD100 sens_D4_LDnone rice_pruned rice_pruned2 rice_pruned3"

cd "$SENS_DIR"

for ds in $DATASETS; do
  echo "=== Processing $ds ==="
  
  # Skip if no BED file
  [ -f "${ds}.bed" ] || { echo "  No BED file, skipping"; continue; }
  
  # PCA (2 components) - skip if already done
  if [ ! -f "${ds}_pca.eigenvec" ]; then
    echo "  Running PCA..."
    plink1.9 --bfile "$ds" --pca 2 --out "${ds}_pca" 2>&1 | tail -1
  else
    echo "  PCA already done"
  fi
  
  # Allele frequencies
  if [ ! -f "${ds}_freq.frq" ]; then
    echo "  Running freq..."
    plink1.9 --bfile "$ds" --freq --out "${ds}_freq" 2>&1 | tail -1
  else
    echo "  Freq already done"
  fi
  
  # HWE
  if [ ! -f "${ds}_hwe.hwe" ]; then
    echo "  Running HWE..."
    plink1.9 --bfile "$ds" --hardy --out "${ds}_hwe" 2>&1 | tail -1
  else
    echo "  HWE already done"
  fi
done

echo ""
echo "=== Computing LD decay for LD dimension ==="
# LD decay: only for LD dimension + LDnone
LD_DS="sens_D4_LDnone sens_D4_LD08 sens_D4_LD05 sens_D4_LD100"

for ds in $LD_DS; do
  echo "=== LD: $ds ==="
  [ -f "${ds}.bed" ] || continue
  
  # Get chromosome list from BIM
  CHROMS=$(awk '{print $1}' "${ds}.bim" | sort -u)
  N_SNPS=$(wc -l < "${ds}.bim")
  echo "  SNPs: $N_SNPS"
  
  # For each chromosome, compute LD within 1000kb windows
  # Use --r2 with --ld-window-kb to limit distance
  # Output all pairs (ld-window-r2 0 means no threshold)
  if [ ! -f "${ds}_ld.ld" ]; then
    echo "  Running LD (may take time for large datasets)..."
    # Use a 1000kb window to limit comparisons
    plink1.9 --bfile "$ds" \
      --r2 \
      --ld-window 99999 \
      --ld-window-kb 1000 \
      --ld-window-r2 0 \
      --out "${ds}_ld" 2>&1 | tail -1
    echo "  LD file size: $(wc -c < ${ds}_ld.ld) bytes, lines: $(wc -l < ${ds}_ld.ld)"
  else
    echo "  LD already done"
    echo "  LD file size: $(wc -c < ${ds}_ld.ld) bytes, lines: $(wc -l < ${ds}_ld.ld)"
  fi
done

echo ""
echo "=== ALL DONE ==="
ls -la "$SENS_DIR"/*_pca.eigenvec "$SENS_DIR"/*_freq.frq "$SENS_DIR"/*_hwe.hwe "$SENS_DIR"/*_ld.ld 2>/dev/null
