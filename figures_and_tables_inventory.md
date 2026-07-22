# Supplementary Materials Inventory
## Manuscript: Evaluating the Impact of Quality-Control Filtering and Population Structure Inference on Population-Genetic Analysis

---

## Main Figures

| Figure | File | Title |
|--------|------|-------|
| 1 | `Figure_1_workflow.pdf` | Reproducible population-genomic analysis framework. Schematic overview of the stepwise QC sensitivity design, population-genetic analyses, sensitivity comparison, and open-source workflow. |
| 2 | `Figure_6_QC_impact.pdf` | Impact of QC filtering on dataset composition. SNP retention, sample retention, and missingness summary across all 16 QC scenarios organized by filtering dimension. |
| 3 | `Figure_8_confusion_matrix.pdf` | Cluster assignment consistency across QC scenarios. 5 × 5 confusion matrices comparing ADMIXTURE K=5 assignments for each scenario against the baseline pipeline, with percent agreement annotated. |
| 4 | `Figure_9_ARI_NMI.pdf` | Pairwise clustering agreement across all 17 datasets. (Left) Adjusted Rand Index (ARI) and (Right) Normalized Mutual Information (NMI) heatmaps. Values range from 0 (no agreement) to 1 (identical clustering after label alignment). |
| 5 | `Figure_10_PCA_correlations.pdf` | PCA similarity across QC scenarios. (Left) PC1 and (Right) PC2 Pearson correlation matrices between all scenario pairs, computed on shared samples. |
| 6 | `Figure_11_genetic_diversity.pdf` | Genetic diversity metrics across QC scenarios. Mean PIC, expected heterozygosity (He), observed heterozygosity (Ho), Shannon diversity index, MAF, and SNP count for each scenario, faceted by filtering dimension. Red dashed line indicates baseline value. |
| 7 | `Figure_12_FST_comparison.pdf` | Effects of QC filtering on mean FST. (Left) Barplot of mean FST for all 16 scenarios, colored by filtering dimension. (Right) Line plot showing FST relative to the baseline (dashed red line). |
| 8 | `Figure_13_AMOVA.pdf` | Among-population variance component (ΦST) from AMOVA across QC scenarios. Dashed line indicates baseline value. |
| 9 | `Figure_14_LD_decay.pdf` | Linkage disequilibrium decay curves for five LD pruning strategies. Mean r² binned by physical distance (log scale). Inset table shows number of SNPs retained per scenario. |
| 10 | `Figure_S5_cluster_flow.pdf` | Redistribution of baseline cluster assignments across QC scenarios. Stacked barplots show the percentage of each baseline cluster assigned to each new cluster under each scenario. |
| 11 | `Figure_16_recommendations.pdf` | Best-practice recommendations dashboard. SNP retention, mean FST, ARI vs. baseline, and short-range LD (0–10 kb mean r²) across all QC thresholds. Dashed red lines indicate baseline or recommended threshold. |

---

## Supplementary Figures

| Figure | File | Title |
|--------|------|-------|
| S1 | `Figure_PCA_structure_validation.pdf` | Principal component analysis of the Rice Diversity Panel 1 (379 accessions, baseline pipeline). PC1 vs. PC2 colored by ADMIXTURE K=5 cluster assignment, showing separation of the five major rice subpopulations. |
| S2 | `Figure_S2_sensitivity_CV_curves.pdf` | ADMIXTURE cross-validation error curves for K = 1–10 across all 16 QC scenarios. CV error decreases rapidly to K = 2 and plateaus around K = 5–7. |
| S3 | `Figure_S3_sensitivity_PCA_scree.pdf` | PCA scree plots (top 10 eigenvalues) for all 16 scenarios. The first two PCs capture the majority of variance in all scenarios. |
| S4 | `Figure_S4_sensitivity_summary.pdf` | Integrated sensitivity summary. Combined visualization of sample size, SNP count, FST, CV error at K=5, and Procrustes similarity across all 16 scenarios. |
| S5 | `Figure_7_cluster_stability.pdf` | Per-sample cluster assignment agreement with the baseline. Distribution of agreement scores (proportion of scenarios assigning each accession to the same cluster as the baseline). |

---

## Supplementary Tables

| Table | File | Title |
|-------|------|-------|
| S1 | `Table_S1_sensitivity_counts.csv` | Dataset composition after each QC scenario. Number of accessions, SNPs, retained SNPs after LD pruning, and SNP missingness rate per scenario. |
| S2 | `Table_S2_sensitivity_CV_K5.csv` | ADMIXTURE cross-validation error at K = 5 for all 16 scenarios. |
| S3 | `Table_S3_PCA_procrustes.csv` | Procrustes similarity between PCA loadings of each scenario and the baseline. |
| S4 | `Table_S4_sensitivity_master.csv` | Master table combining sample size, SNP count, FST, CV error, diversity metrics, and cluster agreement for all 16 scenarios. |
| S5 | `Table_S5_cluster_stability.csv` | Per-accession cluster assignments (after label alignment) for all 17 datasets and the overall agreement score. |
| S6 | `Table_S6_cluster_flow.csv` | Cluster flow table. Percentage of each baseline cluster reassigned to each cluster under each QC scenario, with absolute counts. |
| S7 | `Table_S7_pairwise_ARI.csv` | Pairwise Adjusted Rand Index (ARI) matrix between all 17 datasets. |
| S8 | `Table_S8_pairwise_NMI.csv` | Pairwise Normalized Mutual Information (NMI) matrix between all 17 datasets. |
| S9 | `Table_S9_PCA_correlations.csv` | Pearson correlation of PC1 and PC2 loadings between all scenario pairs. |
| S10 | `Table_S10_genetic_diversity.csv` | Mean PIC, expected heterozygosity (He), observed heterozygosity (Ho), Shannon diversity index, and MAF for all 16 scenarios. |
| S11 | `Table_S11_FST_comparison.csv` | Mean FST, standard deviation, sample size, and SNP count for all 16 scenarios. |
| S12 | `Table_S12_AMOVA.csv` | AMOVA results: among-population and within-population variance components, ΦST, and percentage of variance explained by population structure for all 16 scenarios. |
| S13 | `Table_S13_LD_decay.csv` | LD decay data. Mean r² per distance bin (0–10, 10–50, 50–100, 100–200, 200–500, 500–1,000, 1,000–5,000 kb) and number of SNP pairs per bin for five LD pruning strategies. |
| S14 | `Table_S14_excluded_accessions.csv` | Excluded accessions. Complete list of 34 accessions removed during QC, with NSFTV ID, accession name, country, subpopulation, and exclusion step. |
| S15 | `Table_S15_recommendations.csv` | Best-practice recommendations. Recommended thresholds for sample missingness, MAF, SNP missingness, and LD pruning with empirical impact on SNP retention, FST, cluster stability, and LD decay. |

---

## File Inventory Summary

| Type | Count |
|------|-------|
| Main Figures | 11 |
| Supplementary Figures | 5 |
| Supplementary Tables | 15 |
| **Total** | **31** |
