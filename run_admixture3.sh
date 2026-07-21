#!/bin/bash
cd /root
/opt/miniconda/bin/admixture --cv rice_pruned3.bed 1 | tee log3_1.out
/opt/miniconda/bin/admixture --cv rice_pruned3.bed 2 | tee log3_2.out
/opt/miniconda/bin/admixture --cv rice_pruned3.bed 3 | tee log3_3.out
/opt/miniconda/bin/admixture --cv rice_pruned3.bed 4 | tee log3_4.out
/opt/miniconda/bin/admixture --cv rice_pruned3.bed 5 | tee log3_5.out
/opt/miniconda/bin/admixture --cv rice_pruned3.bed 6 | tee log3_6.out
/opt/miniconda/bin/admixture --cv rice_pruned3.bed 7 | tee log3_7.out
/opt/miniconda/bin/admixture --cv rice_pruned3.bed 8 | tee log3_8.out
/opt/miniconda/bin/admixture --cv rice_pruned3.bed 9 | tee log3_9.out
/opt/miniconda/bin/admixture --cv rice_pruned3.bed 10 | tee log3_10.out
echo "ALL DONE"