#!/bin/bash
cd /home/yehia/hdra/HDRA-G6-4-RDP1-RDP2-NIAS

echo "=== RDP1 samples ==="
awk -F"\t" '$5 == "RDP1"' HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv | wc -l

echo "=== Subpop distribution ==="
awk -F"\t" '$5 == "RDP1" {print $6}' HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv | sort | uniq -c | sort -rn

echo "=== Overlap with 44K ==="
awk -F"\t" '$5 == "RDP1" {print $3}' HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv > /tmp/hdra_rdp1.txt
awk '{print "NSFTV_"$1}' /home/yehia/sensitivity/rice_pruned3.fam > /tmp/44k_ids.txt
comm -12 <(sort /tmp/hdra_rdp1.txt) <(sort /tmp/44k_ids.txt) | wc -l

echo "=== Total HDRA samples ==="
wc -l HDRA-G6-4-RDP1-RDP2-NIAS.fam