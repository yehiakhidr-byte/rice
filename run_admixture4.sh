#!/bin/bash
# ADMIXTURE K=1-12 with 10 replicates each
# Using rice_pruned3 (379 samples x 1187 LD-pruned SNPs)
set -e

cd /root
DATA=rice_pruned3
MIN_K=1
MAX_K=12
NREP=10

rm -rf admixture_results
mkdir -p admixture_results

for K in $(seq $MIN_K $MAX_K); do
    mkdir -p admixture_results/K$K
    best_ll=-1e15
    best_rep=""
    best_cv=""

    for rep in $(seq 1 $NREP); do
        seed=$((rep * 10000 + K))
        logfile="admixture_results/K${K}/log_K${K}_rep${rep}.out"

        /usr/local/bin/admixture --cv --seed=$seed ${DATA}.bed $K > $logfile 2>&1

        # Rename output files
        mv ${DATA}.${K}.Q  admixture_results/K${K}/K${K}_rep${rep}.Q 2>/dev/null || true
        mv ${DATA}.${K}.P  admixture_results/K${K}/K${K}_rep${rep}.P 2>/dev/null || true

        # Extract loglikelihood from log file
        ll=$(grep "Loglikelihood:" $logfile | tail -1 | awk '{print $2}')
        cv=$(grep "CV error" $logfile | awk '{print $4}')

        echo "K=$K rep=$rep seed=$seed LL=$ll CV=$cv"

        # Track best replicate (highest loglikelihood) - use awk for float compare
        better=$(awk "BEGIN {print ($ll > $best_ll) ? 1 : 0}")
        if [ "$better" = "1" ]; then
            best_ll=$ll
            best_rep=$rep
            best_cv=$cv
        fi
    done

    echo "K=$K best rep=$best_rep LL=$best_ll CV=$best_cv"
    echo "$best_ll $best_cv $best_rep" > admixture_results/K${K}/best.txt

    # Copy best Q and P to top level
    cp admixture_results/K${K}/K${K}_rep${best_rep}.Q admixture_results/rice_pruned3.${K}.Q
    cp admixture_results/K${K}/K${K}_rep${best_rep}.P admixture_results/rice_pruned3.${K}.P
done

# Summary of all best CV errors
echo "=== Final CV Error Summary ==="
for K in $(seq $MIN_K $MAX_K); do
    cv=$(grep "CV error" admixture_results/K${K}/log_K${K}_rep*.out | head -1 | awk '{print $4}')
    echo "K=$K CV=$cv"
done
