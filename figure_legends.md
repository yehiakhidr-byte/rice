# Figure and Table Legends

---

## Main Figures

### Figure 1. Experimental workflow for evaluating QC sensitivity in the Rice Diversity Panel 1.
The stepwise framework used in this study. Starting from the raw RDP1 dataset (413 accessions, 26,474 SNPs), four QC dimensions were independently varied: sample missingness (--mind), minor allele frequency (--maf), SNP missingness (--geno), and LD pruning (--indep-pairwise). A total of 16 scenarios were generated (5 + 4 + 4 + 5, each with the baseline in one dimension). For each scenario, population-genetic analyses (PCA, ADMIXTURE K = 1–10, FST, AMOVA, genetic diversity, LD decay) were performed. Sensitivity was assessed using confusion matrices, Adjusted Rand Index (ARI), Normalized Mutual Information (NMI), and PC correlations. Results were synthesized into evidence-based recommendations. All code, parameter files, and documentation are publicly available.

### Figure 2. Impact of quality-control filtering on dataset composition.
Barplots showing the number of SNPs retained (top) and number of accessions retained (bottom) for each of the 16 QC scenarios, grouped by filtering dimension (blue: sample filtering; red: MAF; green: SNP missingness; purple: LD pruning). The baseline scenario (mind 0.10, MAF 0.05, geno 0.05, LD 50/5/0.2) is highlighted with a dashed red line.

### Figure 3. Cluster assignment consistency across QC scenarios.
5 × 5 confusion matrices comparing ADMIXTURE K = 5 cluster assignments for each non-baseline scenario (rows: baseline clusters C1–C5; columns: scenario clusters). Cell values show the percentage of accessions from each baseline cluster assigned to each scenario cluster after label alignment. The overall agreement rate is annotated above each matrix. Colors range from white (0%) to dark blue (100%). Misassignments predominantly involve adjacent clusters (e.g., C1 ↔ C5), consistent with the genetic continuum among rice subpopulations.

### Figure 4. Effects of quality-control filtering on population differentiation (FST).
Barplot of mean pairwise FST (Weir & Cockerham, 1984) for all 16 QC scenarios, colored by filtering dimension (blue: sample; red: MAF; green: SNP missingness; purple: LD pruning). The dashed red line marks the baseline FST (0.423). FST values range from 0.293 (mind 0.02) to 0.609 (no LD pruning). The stability of FST across MAF, GENO, and moderate LD scenarios (0.39–0.45) indicates that population-differentiation estimates are robust to reasonable QC choices. The inflation to 0.609 in the absence of LD pruning demonstrates that correlated markers systematically bias FST upward.

### Figure 5. Linkage disequilibrium decay under different LD pruning strategies.
Mean r² between SNP pairs plotted against physical distance (kb, log scale) for five LD pruning scenarios: no pruning (26,474 SNPs), r² < 0.8 (9,609 SNPs), r² < 0.5 (4,278 SNPs), 50/5/0.2 — baseline (1,187 SNPs), and 100/10/0.2 (822 SNPs). Vertical axis shows mean r² per distance bin (bins: 0–10, 10–50, 50–100, 100–200, 200–500, 500–1,000, 1,000–5,000 kb). Without pruning, mean r² at 0–10 kb is 0.43, with 48.1% of SNP pairs exceeding r² > 0.2. The baseline window reduces mean r² to 0.074 at 0–10 kb (5.8% > 0.2). The qualitative pattern — rapid decay within 100 kb followed by gradual approach to background — is consistent with the original RDP1 report (Zhao et al., 2011).

### Figure 6. Best-practice recommendations for reproducible population-genomic analysis of RDP1.
Four-panel dashboard summarizing the impact of each QC threshold on key analytical endpoints. **(a)** SNP retention (number of SNPs). **(b)** Mean FST (dashed red line: baseline). **(c)** Cluster stability (ARI vs. baseline; dashed line: ARI = 0.9). **(d)** Short-range LD (mean r² at 0–10 kb; dashed line: baseline). Colors indicate filtering dimension. Recommended thresholds: mind 0.10, MAF 0.05, geno 0.05, LD pruning 50/5/0.2.

---

## Supplementary Figures

### Figure S1. Principal component analysis of the Rice Diversity Panel 1 (baseline pipeline, 379 accessions, 1,187 SNPs).
PC1 versus PC2, with each accession colored by its ADMIXTURE K = 5 cluster assignment after label alignment. The five major rice subpopulations (indica, temperate japonica, tropical japonica, aus, aromatic) are separated along the first two principal components, consistent with the established population structure of RDP1 (Zhao et al., 2011). Admixed accessions occupy intermediate positions. Percent variance explained: PC1 = 28.3%, PC2 = 10.1%.

### Figure S2. ADMIXTURE cross-validation error curves for K = 1–10 across all 16 QC scenarios.
Each panel shows CV error (y-axis) as a function of K (x-axis) for one scenario. CV error decreases rapidly from K = 1 to K = 2, plateaus around K = 5–7, and increases slightly at K > 7. The CV minimum occurs at K = 5 for 10 of 16 scenarios and at K = 4 for the remaining 6. The qualitative shape is consistent across all QC scenarios, indicating that the optimal K is robust to filtering decisions.

### Figure S3. Pairwise clustering agreement across all 17 datasets.
Heatmaps of the Adjusted Rand Index (left) and Normalized Mutual Information (right) between all pairs of datasets (16 QC scenarios + baseline). Values range from 0 (no agreement) to 1 (identical clustering after label alignment). Within-dimension agreement is highest for MAF and GENO scenarios (ARI > 0.95), intermediate for LD scenarios, and lowest for sample-filtering scenarios (ARI range: 0.10–0.58).

### Figure S4. PCA similarity across QC scenarios.
Pearson correlation matrices for PC1 (left) and PC2 (right) loadings between all scenario pairs, computed on the intersection of shared samples. PC1 correlations exceed 0.99 for all MAF, GENO, and LD scenarios versus the baseline. The only notable deviation is mind = 0.02 (PC1 r = 0.82). PC2 correlations are slightly more variable but remain above 0.90.

### Figure S5. Distribution of per-accession cluster assignment agreement with the baseline.
Histogram showing, for each of the 379 accessions, the proportion of the 15 non-baseline scenarios in which its cluster assignment (after label alignment) matches the baseline. Mean agreement = 93.1%. Of all accessions, 43.3% are fully stable (agree = 1.0) and 81% have agreement ≥ 0.9. Only 12.4% of accessions have agreement < 0.8, predominantly those affected by the most stringent sample-filtering scenarios (mind 0.02, mind 0.05).

---

## Main Tables

### Table 1. Accessions excluded during quality-control filtering.
Complete list of 34 accessions removed from the original 413-accession RDP1 panel during the baseline QC pipeline, with NSFTV identifier, accession name, country of origin, subpopulation classification, and exclusion step. Four accessions were removed by --mind 0.10 (>10% genotype missingness); 30 were lost during the stepwise SNP filtering pipeline (MAF → GENO → LD) when previously undetected sample-level missingness exceeded the implied threshold. The excluded accessions are distributed roughly proportionally across subpopulations (ADMIX: −8, AROMATIC: −1, AUS: −8, IND: −10, TEJ: −3, TRJ: −4; overall exclusion rate = 8.2%). Country representation is also preserved (top 8 countries unchanged).

### Table 2. Best-practice recommendations for population-genomic analysis of the Rice Diversity Panel 1.
Recommended QC thresholds with empirical impact on five key analytical endpoints: SNP retention, sample retention, genetic diversity (PIC), population differentiation (FST), and cluster stability (ARI vs. baseline). Values are shown for each recommended threshold and for the nearest alternatives to illustrate the sensitivity range. The recommended pipeline is: --mind 0.10 → --maf 0.05 → --geno 0.05 → --indep-pairwise 50 5 0.2, yielding 379 accessions × 1,187 SNPs with stable population-genetic estimates.

---

## Supplementary Tables

### Table S1. Dataset composition after each QC scenario.
Number of accessions and SNPs retained under each of the 16 QC scenarios. Columns: Scenario ID, QC dimension, parameter value, number of accessions, number of SNPs before LD pruning, number of SNPs after LD pruning (if applicable), and SNP missingness rate.

### Table S2. ADMIXTURE cross-validation error at K = 5 for all 16 scenarios.

### Table S3. Procrustes similarity between PCA loadings.
Comparison of PCA loadings from each scenario against the baseline, measured by Procrustes sum of squares (lower = more similar).

### Table S4. Master sensitivity table.
Combined results for all 16 scenarios: sample size, SNP count, mean FST, ADMIXTURE CV error at K = 5, mean PIC, mean He, mean Ho, Shannon diversity index, and cluster agreement with baseline.

### Table S5. Per-accession cluster assignments for all 17 datasets.
For each of the 379 accessions (rows), the ADMIXTURE K = 5 cluster assignment (after label alignment) in the baseline and all 15 non-baseline scenarios, plus the overall agreement score.

### Table S6. Cluster flow between baseline and each scenario.
For each baseline cluster and each non-baseline scenario, the number and percentage of accessions reassigned to each new cluster.

### Table S7. Pairwise Adjusted Rand Index (ARI) between all 17 datasets.

### Table S8. Pairwise Normalized Mutual Information (NMI) between all 17 datasets.

### Table S9. PC1 and PC2 Pearson correlation coefficients between all scenario pairs.

### Table S10. Genetic diversity metrics for all 16 scenarios.
Mean PIC, expected heterozygosity (He), observed heterozygosity (Ho), Shannon diversity index, mean MAF, and number of SNPs per scenario.

### Table S11. Mean FST and AMOVA results for all 16 scenarios.
Mean pairwise FST, FST standard deviation, AMOVA among-population variance (ΦST), within-population variance, and percentage of variance explained by population structure.

### Table S12. LD decay data for five LD pruning strategies.
Mean r² per physical distance bin (0–10, 10–50, 50–100, 100–200, 200–500, 500–1,000, 1,000–5,000 kb) and number of SNP pairs per bin.

### Table S13. Complete excluded accession list with passport data.
(Identical to Table 1, provided as a standalone CSV for computational access.)

### Table S14. Best-practice recommendation details.
(Identical to Table 2, provided as a standalone CSV for computational access.)