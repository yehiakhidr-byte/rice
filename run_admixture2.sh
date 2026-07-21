#!/bin/bash
cd /root
/opt/miniconda/bin/admixture --cv rice_pruned2.bed 1 | tee log2_1.out
/opt/miniconda/bin/admixture --cv rice_pruned2.bed 2 | tee log2_2.out
/opt/miniconda/bin/admixture --cv rice_pruned2.bed 3 | tee log2_3.out
/opt/miniconda/bin/admixture --cv rice_pruned2.bed 4 | tee log2_4.out
/opt/miniconda/bin/admixture --cv rice_pruned2.bed 5 | tee log2_5.out
/opt/miniconda/bin/admixture --cv rice_pruned2.bed 6 | tee log2_6.out
/opt/miniconda/bin/admixture --cv rice_pruned2.bed 7 | tee log2_7.out
/opt/miniconda/bin/admixture --cv rice_pruned2.bed 8 | tee log2_8.out
/opt/miniconda/bin/admixture --cv rice_pruned2.bed 9 | tee log2_9.out
/opt/miniconda/bin/admixture --cv rice_pruned2.bed 10 | tee log2_10.out
echo "ALL DONE"