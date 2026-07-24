#!/bin/bash
# HDRA RDP1 Sensitivity Pipeline
set -e
cd /home/yehia/hdra

HDRA_DIR="/home/yehia/hdra/HDRA-G6-4-RDP1-RDP2-NIAS"
OUTDIR="/home/yehia/hdra_rdp1"
mkdir -p $OUTDIR

echo "=== Step 1: Extract RDP1 samples from HDRA ==="

# Create keep list from sample_map (column 5 = RDP1, column 1 = FAM ID)
awk -F"\t" '$5 == "RDP1" {print $1"\t"$1}' $HDRA_DIR/HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv > /tmp/rdp1_keep.txt
echo "RDP1 samples to extract: $(wc -l < /tmp/rdp1_keep.txt)"

# Also create NSFTV mapping for later comparison with 44K
awk -F"\t" '$5 == "RDP1" {print $1"\t"$3}' $HDRA_DIR/HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv > /tmp/rdp1_nsftv_map.txt
echo "NSFTV mapping created: $(wc -l < /tmp/rdp1_nsftv_map.txt)"

# Extract with PLINK
plink1.9 --bfile $HDRA_DIR/HDRA-G6-4-RDP1-RDP2-NIAS \
  --keep /tmp/rdp1_keep.txt \
  --make-bed \
  --out $OUTDIR/hdra_rdp1_raw \
  2>&1 | tail -5

echo ""
echo "=== Step 2: QC on raw HDRA RDP1 ==="
# Check missingness
plink1.9 --bfile $OUTDIR/hdra_rdp1_raw \
  --missing \
  --out $OUTDIR/hdra_rdp1_raw 2>&1 | tail -3

# Sample missingness summary
echo "Sample missingness (first 10):"
head -11 $OUTDIR/hdra_rdp1_raw.imiss

echo ""
echo "=== Raw dataset stats ==="
echo "Samples: $(wc -l < $OUTDIR/hdra_rdp1_raw.fam)"
echo "SNPs: $(wc -l < $OUTDIR/hdra_rdp1_raw.bim)"