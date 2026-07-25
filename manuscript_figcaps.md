## Figure Captions

**Figure 1.** Workflow diagram of the complete cross-platform QC sensitivity analysis pipeline. Input data (44K and HDRA) follow identical stepwise filtering protocols with parallel sensitivity scenarios, ADMIXTURE, PCA, post-processing, and cross-platform comparison stages.

**Figure 2.** Effects of QC filtering on SNP and sample retention. (a) 44K platform: sample retention by mind threshold and SNP retention by LD pruning threshold. (b) HDRA platform: SNP retention by MAF and GENO thresholds. The recommended baseline thresholds are highlighted.

**Figure 3.** Cross-platform comparison of QC sensitivity. (a) FST values across all scenarios for 44K (blue) and HDRA (red). The MAF dimension shows the largest divergence between platforms. (b) Cluster assignment agreement with the baseline scenario. HDRA shows systematically lower agreement, especially in the sample and MAF dimensions. (c) Genetic diversity (He) across scenarios. HDRA He spans a 6.3-fold range driven by MAF filtering, while 44K He is nearly constant.

**Figure 4.** PCA robustness across platforms. PC1 correlation between each scenario and the baseline is plotted for 44K and HDRA. Both platforms maintain PC1 ρ ≥ 0.97 across all scenarios, indicating that the primary axis of genetic variation is insensitive to QC platform.

**Figure 5.** FST as a function of MAF threshold in HDRA. FST increases from 0.076 (no MAF filter) to 0.443 (MAF = 0.10), demonstrating that MAF filtering is the dominant determinant of differentiation estimates in high-density data.

**Figure 6.** LD-driven FST inflation. The ratio of unpruned-to-baseline FST is 1.44× in 44K and 1.81× in HDRA, reflecting the greater residual LD in higher-density data after equivalent pruning parameters.

**Figure 7.** ADMIXTURE analysis of 379 44K accessions based on 1,187 LD-pruned SNPs. (a) Cross-validation error as a function of K (K = 1–12). The red point highlights K = 5. (b) Ancestry proportions at K = 2. (c) Ancestry proportions at K = 5.

**Figure 8.** DAPC analysis of 379 44K accessions. (a) Scatter plot at K = 2. (b) Scatter plot at K = 5 showing five genetic clusters. (c) Posterior membership probabilities at K = 5.

**Figure 9.** Kinship matrix heatmap of 379 44K accessions based on 26,474 QC-passing SNPs. The centered IBS kinship matrix is visualized with hierarchical clustering (Euclidean distance, complete linkage).

**Figure 10.** LD decay curves on chromosome 1 for both platforms under different pruning strategies. (a) 44K: progressive pruning reduces short-range LD from r² = 0.43 (unpruned) to 0.074 (50/5/0.2). (b) HDRA: comparable pruning effects with higher residual LD in the unpruned dataset.
