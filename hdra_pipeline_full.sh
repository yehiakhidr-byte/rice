#!/bin/bash
# HDRA RDP1 Sensitivity Pipeline - Full
# Generates all 16 scenarios, runs PCA, ADMIXTURE
set -e
cd /home/yehia/hdra_rdp1

BASE="hdra_rdp1_raw"

echo "============================================"
echo "HDRA RDP1 Sensitivity Pipeline"
echo "406 samples x 700K SNPs"
echo "============================================"

# =========================================================
# SCENARIO DIMENSIONS (matching the 44K analysis)
# =========================================================
# D1: Sample filtering (mind)
# S0: none, S1: 0.10_after, S2: 0.05, S3: 0.10 (baseline), S4: 0.02
# D2: MAF
# M0: none, M1: 0.01, M2: 0.03, M3: 0.05 (baseline), M4: 0.10
# D3: GENO
# G0: none, G1: 0.20, G2: 0.10, G3: 0.05 (baseline), G4: 0.02
# D4: LD pruning
# L0: none, L1: 0.8, L2: 0.5, L3: 0.2 (baseline), L4: 100/10/0.2

run_pipeline() {
  local tag=$1
  local mind=$2
  local maf=$3
  local geno=$4
  local ld_window=$5
  local ld_step=$6
  local ld_r2=$7
  
  echo ""
  echo "=== Processing scenario: $tag ==="
  
  # Check if complete
  if [ -f "${tag}_K5.Q" ] && [ -f "${tag}_pca.eigenvec" ]; then
    echo "  Already complete, skipping"
    return
  fi
  
  # Step 1: Sample filtering
  if [ "$mind" = "none" ]; then
    cp ${BASE}.bed ${tag}.bed
    cp ${BASE}.bim ${tag}.bim
    cp ${BASE}.fam ${tag}.fam
  else
    echo "  mind = $mind"
    plink1.9 --bfile ${BASE} --mind $mind --make-bed --out ${tag} 2>&1 | tail -1
  fi
  
  # Step 2: MAF filtering  
  if [ "$maf" != "none" ]; then
    echo "  maf = $maf"
    plink1.9 --bfile ${tag} --maf $maf --make-bed --out ${tag}_tmp1 2>&1 | tail -1
    mv ${tag}_tmp1.bed ${tag}.bed
    mv ${tag}_tmp1.bim ${tag}.bim
    mv ${tag}_tmp1.fam ${tag}.fam
  fi
  
  # Step 3: GENO filtering
  if [ "$geno" != "none" ]; then
    echo "  geno = $geno"
    plink1.9 --bfile ${tag} --geno $geno --make-bed --out ${tag}_tmp2 2>&1 | tail -1
    mv ${tag}_tmp2.bed ${tag}.bed
    mv ${tag}_tmp2.bim ${tag}.bim
    mv ${tag}_tmp2.fam ${tag}.fam
  fi
  
  # Step 4: LD pruning
  if [ "$ld_r2" != "none" ]; then
    echo "  LD pruning: window=$ld_window step=$ld_step r2=$ld_r2"
    plink1.9 --bfile ${tag} --indep-pairwise $ld_window $ld_step $ld_r2 \
      --out ${tag}_ld 2>&1 | tail -1
    n_pruned=$(wc -l < ${tag}_ld.prune.in 2>/dev/null || echo 0)
    echo "  SNPs after LD pruning: $n_pruned"
    plink1.9 --bfile ${tag} --extract ${tag}_ld.prune.in --make-bed \
      --out ${tag} 2>&1 | tail -1
  fi
  
  # Log dataset stats
  n_snps=$(wc -l < ${tag}.bim)
  n_samps=$(wc -l < ${tag}.fam)
  echo "  Dataset: $n_samps samples x $n_snps SNPs"
  echo "$tag,$n_samps,$n_snps" >> hdra_scenario_counts.csv
  
  # Step 5: PCA
  echo "  Running PCA..."
  plink1.9 --bfile ${tag} --pca 10 --out ${tag}_pca 2>&1 | tail -2
  
  # Step 6: Allele frequency
  plink1.9 --bfile ${tag} --freq --out ${tag}_freq 2>&1 | tail -1
  
  # Step 7: HWE
  plink1.9 --bfile ${tag} --hardy --out ${tag}_hwe 2>&1 | tail -1
  
  # Step 8: ADMIXTURE K=1-10
  echo "  Running ADMIXTURE K=1-10..."
  for K in 1 2 3 4 5 6 7 8 9 10; do
    if [ ! -f "${tag}_K${K}.Q" ]; then
      echo "    K=$K..."
      admixture --cv ${tag}.bed $K -j4 2>&1 | tail -1
      mv ${tag}.${K}.Q ${tag}_K${K}.Q
      mv ${tag}.${K}.P ${tag}_K${K}.P
    else
      echo "    K=$K exists, skipping"
    fi
  done
  
  # Extract CV error for K=5
  grep "CV" ${tag}_K5.log 2>/dev/null | awk '{print $NF}' > ${tag}_cv5.txt
  
  echo "  Scenario $tag complete"
}

# =========================================================
# Generate all 16 scenarios
# =========================================================
mkdir -p logs

echo "Creating baseline SNP count file..."
echo "scenario,samples,snps" > hdra_scenario_counts.csv

# D1: Sample filtering dimension
# Start from raw for each
for mind_val in "none" "0.1" "0.05" "0.02"; do
  tag="hdra_S_${mind_val}"
  if [ "$mind_val" = "none" ]; then tag="hdra_S_none"; fi
  if [ "$mind_val" = "0.1" ]; then tag="hdra_S_010"; fi
  if [ "$mind_val" = "0.05" ]; then tag="hdra_S_005"; fi
  if [ "$mind_val" = "0.02" ]; then tag="hdra_S_002"; fi
  
  run_pipeline "$tag" "$mind_val" "0.05" "0.05" "50" "5" "0.2"
done

# Also add S3 (baseline) = mind 0.10 with full stepwise
tag="hdra_S3_baseline"
if [ ! -f "${tag}_K5.Q" ]; then
  run_pipeline "$tag" "0.1" "0.05" "0.05" "50" "5" "0.2"
fi

# D2: MAF dimension (start from baseline-mind, vary MAF)
BASE_MIND="hdra_S_010"
for maf_val in "none" "0.01" "0.03" "0.05" "0.10"; do
  tag="hdra_M_${maf_val}"
  if [ "$maf_val" = "none" ]; then tag="hdra_M_none"; fi
  if [ "$maf_val" = "0.01" ]; then tag="hdra_M_001"; fi
  if [ "$maf_val" = "0.03" ]; then tag="hdra_M_003"; fi
  if [ "$maf_val" = "0.05" ]; then tag="hdra_M_005"; fi
  if [ "$maf_val" = "0.10" ]; then tag="hdra_M_010"; fi
  
  # For MAF scenarios, copy baseline-mind as starting point
  if [ ! -f "${tag}.bed" ]; then
    cp ${BASE_MIND}.bed ${tag}.bed
    cp ${BASE_MIND}.bim ${tag}.bim
    cp ${BASE_MIND}.fam ${tag}.fam
  fi
  run_pipeline "$tag" "0.1" "$maf_val" "0.05" "50" "5" "0.2"
done

# D3: GENO dimension (start from MAF=0.05 baseline)
for geno_val in "none" "0.20" "0.10" "0.05" "0.02"; do
  tag="hdra_G_${geno_val}"
  if [ "$geno_val" = "none" ]; then tag="hdra_G_none"; fi
  if [ "$geno_val" = "0.20" ]; then tag="hdra_G_020"; fi
  if [ "$geno_val" = "0.10" ]; then tag="hdra_G_010"; fi
  if [ "$geno_val" = "0.05" ]; then tag="hdra_G_005"; fi
  if [ "$geno_val" = "0.02" ]; then tag="hdra_G_002"; fi
  
  if [ ! -f "${tag}.bed" ]; then
    cp ${BASE_MIND}.bed ${tag}.bed
    cp ${BASE_MIND}.bim ${tag}.bim
    cp ${BASE_MIND}.fam ${tag}.fam
  fi
  run_pipeline "$tag" "0.1" "0.05" "$geno_val" "50" "5" "0.2"
done

# D4: LD pruning dimension (start from GENO=0.05 baseline)
for ld_r2 in "none" "0.8" "0.5" "0.2" "custom"; do
  tag="hdra_L_${ld_r2}"
  if [ "$ld_r2" = "none" ]; then tag="hdra_L_none"; fi
  if [ "$ld_r2" = "0.8" ]; then tag="hdra_L_08"; fi
  if [ "$ld_r2" = "0.5" ]; then tag="hdra_L_05"; fi
  if [ "$ld_r2" = "0.2" ]; then tag="hdra_L_02"; fi
  if [ "$ld_r2" = "custom" ]; then tag="hdra_L_100"; fi
  
  if [ ! -f "${tag}.bed" ]; then
    cp ${BASE_MIND}.bed ${tag}.bed
    cp ${BASE_MIND}.bim ${tag}.bim
    cp ${BASE_MIND}.fam ${tag}.fam
  fi
  
  if [ "$ld_r2" = "custom" ]; then
    run_pipeline "$tag" "0.1" "0.05" "0.05" "100" "10" "0.2"
  else
    run_pipeline "$tag" "0.1" "0.05" "0.05" "50" "5" "$ld_r2"
  fi
done

echo ""
echo "============================================"
echo "All 16 scenarios generated"
cat hdra_scenario_counts.csv
echo "============================================"
echo "Pipeline complete - next: FST, AMOVA, diversity, LD decay"