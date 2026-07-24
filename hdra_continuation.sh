#!/bin/bash
# HDRA RDP1 Sensitivity Pipeline - Continuation
# Skips degenerate scenarios (e.g., S_002 with only 1 sample)
set -e
cd /home/yehia/hdra_rdp1
BASE="hdra_rdp1_raw"
N_CORES=4

gen_tag() {
  local d=$1 v=$2
  case $d in
    S) case $v in none) echo hdra_S_none;; 0.1) echo hdra_S_010;; 0.05) echo hdra_S_005;; 0.02) echo hdra_S_002;; esac ;;
    M) case $v in none) echo hdra_M_none;; 0.01) echo hdra_M_001;; 0.03) echo hdra_M_003;; 0.05) echo hdra_M_005;; 0.10) echo hdra_M_010;; esac ;;
    G) case $v in none) echo hdra_G_none;; 0.20) echo hdra_G_020;; 0.10) echo hdra_G_010;; 0.05) echo hdra_G_005;; 0.02) echo hdra_G_002;; esac ;;
    L) case $v in none) echo hdra_L_none;; 0.8) echo hdra_L_08;; 0.5) echo hdra_L_05;; 0.2) echo hdra_L_02;; 100_10_0.2) echo hdra_L_100;; esac ;;
  esac
}

run_one() {
  local t=$1 m=$2 a=$3 g=$4 w=$5 s=$6 r=$7
  echo "[$(date)] Starting: $t"
  
  # Check if already complete
  if [ -f "${t}_K10.Q" ] && [ -f "${t}_pca.eigenvec" ]; then
    echo "  SKIP $t (already complete)"
    echo "$t,skip" >> p.csv
    return
  fi
  
  # Apply filtering
  if [ ! -f "${t}.bed" ]; then
    local cmd="plink1.9 --bfile ${BASE}"
    [ "$m" != none ] && cmd="$cmd --mind $m"
    [ "$a" != none ] && cmd="$cmd --maf $a"
    [ "$g" != none ] && cmd="$cmd --geno $g"
    cmd="$cmd --make-bed --out ${t}"
    eval $cmd 2>&1 | tail -1
  fi
  
  n_samps=$(wc -l < ${t}.fam)
  n_snps=$(wc -l < ${t}.bim)
  echo "  After filtering: $n_samps samples, $n_snps SNPs"
  
  # Skip if too few samples or SNPs
  if [ "$n_samps" -lt 5 ] || [ "$n_snps" -lt 10 ]; then
    echo "  SKIP $t (degenerate: $n_samps samples, $n_snps SNPs)"
    echo "$t,fail($n_samps x $n_snps)" >> p.csv
    return
  fi
  
  # LD prune
  if [ "$r" != none ]; then
    if [ ! -f "${t}.prune.in" ]; then
      plink1.9 --bfile ${t} --indep-pairwise $w $s $r --out ${t} 2>&1 | tail -1
    fi
    if [ ! -f "${t}.prune.in" ]; then
      echo "  SKIP $t (LD pruning failed)"
      echo "$t,fail(LD)" >> p.csv
      return
    fi
    if [ ! -f "${t}.pruned.bed" ]; then
      plink1.9 --bfile ${t} --extract ${t}.prune.in --make-bed --out ${t}.pruned 2>&1 | tail -1
    fi
    b="${t}.pruned"
  else
    b="${t}"
  fi
  
  # Log stats
  local ns=$(wc -l < ${b}.fam)
  local nk=$(wc -l < ${b}.bim)
  echo "  Final: $ns samples x $nk SNPs"
  # Only add to s.csv if not already there
  grep -q "^$t," s.csv 2>/dev/null || echo "$t,$ns,$nk" >> s.csv
  
  # PCA
  [ ! -f "${t}_pca.eigenvec" ] && plink1.9 --bfile $b --pca 10 --out ${t}_pca 2>&1 | tail -1
  
  # Frequencies
  [ ! -f "${t}_freq.frq" ] && plink1.9 --bfile $b --freq --out ${t}_freq 2>&1 | tail -1
  [ ! -f "${t}_hwe.hwe" ] && plink1.9 --bfile $b --hardy --out ${t}_hwe 2>&1 | tail -1
  
  # ADMIXTURE
  for K in 1 2 3 4 5 6 7 8 9 10; do
    if [ ! -f "${t}_K${K}.Q" ]; then
      echo "  ADMIXTURE K=$K..."
      admixture --cv ${b}.bed $K -j${N_CORES} 2>&1 | tail -1
      mv ${b}.${K}.Q ${t}_K${K}.Q 2>/dev/null
      mv ${b}.${K}.P ${t}_K${K}.P 2>/dev/null
    fi
  done
  
  echo "$t,done" >> p.csv
  echo "[$(date)] Done: $t"
}

echo "=== HDRA Pipeline Continuation ==="

# D1: Sample filtering (skip S_002 if degenerate)
for v in none 0.1 0.05; do
  run_one "$(gen_tag S $v)" "$v" 0.05 0.05 50 5 0.2
done

# D2: MAF (start from mind=0.1 baseline)
for v in none 0.01 0.03 0.05 0.10; do
  run_one "$(gen_tag M $v)" 0.1 "$v" 0.05 50 5 0.2
done

# D3: GENO
for v in none 0.20 0.10 0.05 0.02; do
  run_one "$(gen_tag G $v)" 0.1 0.05 "$v" 50 5 0.2
done

# D4: LD
run_one "$(gen_tag L none)" 0.1 0.05 0.05 "" "" none
run_one "$(gen_tag L 0.8)"  0.1 0.05 0.05 50 5 0.8
run_one "$(gen_tag L 0.5)"  0.1 0.05 0.05 50 5 0.5
run_one "$(gen_tag L 0.2)"  0.1 0.05 0.05 50 5 0.2
run_one "$(gen_tag L 100_10_0.2)" 0.1 0.05 0.05 100 10 0.2

echo "=== CONTINUATION COMPLETE ==="
cat s.csv