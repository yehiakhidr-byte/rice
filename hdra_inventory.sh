#!/bin/bash
cd /home/yehia/hdra/HDRA-G6-4-RDP1-RDP2-NIAS

echo "=== RDP1 samples in HDRA ==="
awk -F"\t" '$5 == "RDP1"' HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv | wc -l

echo ""
echo "=== RDP1 FAM IDs (first 10) ==="
awk -F"\t" '$5 == "RDP1" {print $1}' HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv | head -10

echo ""
echo "=== Subpopulation breakdown in HDRA RDP1 ==="
awk -F"\t" '$5 == "RDP1" {print $6}' HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv | sort | uniq -c | sort -rn

echo ""
echo "=== Overlap with 44K RDP1 ==="
awk -F"\t" '$5 == "RDP1" {print $3}' HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv | sort > /tmp/hdra_rdp1_ids.txt
awk '{print "NSFTV_"$1}' /home/yehia/sensitivity/rice_pruned3.fam | sort > /tmp/44k_ids.txt
comm -12 /tmp/hdra_rdp1_ids.txt /tmp/44k_ids.txt | wc -l