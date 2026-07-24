## Materials and Methods

### Model-based ancestry inference (ADMIXTURE)

To complement the STRUCTURE analysis, model-based ancestry estimation was performed using ADMIXTURE v1.3.0 (Alexander et al. 2009) on the LD-pruned dataset of 379 accessions and 1,187 approximately independent SNPs (r² < 0.2 within 50-SNP sliding windows). Analyses were conducted for K = 1 to 12, with 10 independent replicates per K using different random seeds. For each K, the replicate with the highest log-likelihood was retained. Cross-validation (CV) error was estimated using 5-fold CV. The optimal K was assessed by inspecting the CV error curve and comparing results with the STRUCTURE-derived K = 2 and K = 5 models.

### Multivariate clustering (DAPC)

Discriminant Analysis of Principal Components (DAPC; Jombart et al. 2010) was implemented in the adegenet v2.1.11 R package using the same 379 × 1,187 SNP matrix. Missing genotypes were imputed with the mean allele frequency at each locus. Principal components were computed from the genotype matrix, and sequential K-means clustering was performed for up to 20 clusters using the `find.clusters()` function. The optimal number of clusters was selected based on the Bayesian Information Criterion (BIC). To enable direct comparison with the STRUCTURE models, DAPC was also run with K = 2 and K = 5. The first two discriminant functions were visualized.

### Kinship analysis

A centered IBS kinship matrix was computed in PLINK v1.9 (`--make-rel`) using 26,474 QC-passing SNPs across all 379 accessions. The resulting 379 × 379 matrix was visualized as a heatmap with hierarchical clustering (Euclidean distance, complete linkage). The distribution of pairwise kinship coefficients was summarized.

---

## Results

### ADMIXTURE validation of population structure

ADMIXTURE cross-validation errors decreased monotonically from K = 1 (CV = 1.049) to K = 12 (CV = 0.546; Table 5). The largest reduction occurred between K = 1 and K = 2 (ΔCV = −0.238), after which the rate of decline diminished progressively. The CV error at K = 5 was 0.621, with only marginal improvements at higher K values (ΔCV from K = 5 to K = 12 = −0.075; Figure 4a).

At K = 2 (Figure 4b), ADMIXTURE partitioned the 379 accessions into two groups of 146 and 233 individuals, compared to 224 and 155 reported by STRUCTURE at the same K. The differences reflect the use of LD-pruned markers (1,187 SNPs) in ADMIXTURE versus the full marker set (36,901 SNPs) in STRUCTURE. Despite this, the overall pattern of two major genetic groups was consistent between methods.

At K = 5 (Figure 4c), ADMIXTURE ancestry proportions showed close correspondence with STRUCTURE membership coefficients at the same K. A cross-tabulation of dominant ancestry assignments between the two methods showed that 97.6% of accessions were assigned to corresponding subpopulations (Table 6), confirming that the five ancestral populations identified by ADMIXTURE align with the five STRUCTURE-defined subpopulations (aus, indica, temperate japonica, tropical japonica, aromatic). At K = 12 (Supplementary Figure S1), additional subdivisions were detected within the major groups but primarily reflected fine-scale structure.

### DAPC clustering

The optimal number of DAPC clusters based on BIC was five. The five clusters contained 14, 89, 50, 110, and 116 accessions, respectively (Table 7). The first discriminant axis explained 61.8% of the between-group variance, and the second contributed an additional 21.0%, together accounting for 82.8% of the total discrimination.

When DAPC was forced to K = 2 (Figure 5a), the two clusters contained 239 and 140 individuals, showing strong correspondence with ADMIXTURE K = 2: 233 of 233 ADMIXTURE cluster-2 individuals were assigned to DAPC cluster 1, while 140 of 146 ADMIXTURE cluster-1 individuals were assigned to DAPC cluster 2. At K = 5 (Figure 5b), the scatter plot of the first two discriminant functions showed clear separation among all five clusters. DAPC cluster assignments showed 89–100% concordance with ADMIXTURE K = 5 ancestry components after label alignment (Table 6). Posterior membership probabilities exceeded 0.9 for 95% of accessions (Figure 5c), indicating high assignment confidence.

### Kinship relationships

Pairwise kinship coefficients among the 379 accessions ranged from −1.05 to 3.63 (mean = −0.005, SD = 0.809; Table 8). Approximately 48% of pairwise comparisons were negative, indicating pairs less related than expected under random population background. The kinship heatmap (Figure 6) revealed a clear block-diagonal structure consistent with the five genetic clusters identified by STRUCTURE and DAPC, with elevated within-group kinship values and near-zero or negative values between groups. Accessions within the same subpopulation showed substantially higher mean kinship (0.31–0.52) compared to across-subpopulation pairs (mean = −0.12). The proportion of pairwise comparisons exceeding the second-degree relationship threshold (>0.125) was 48.7%, while 39.7% exceeded the first-degree threshold (>0.25), reflecting the presence of closely related accessions within subpopulations.

---

## Table and Figure Captions

**Table 5.** ADMIXTURE cross-validation errors (best of 10 replicates per K) for K = 1–12 based on 379 rice accessions genotyped at 1,187 LD-pruned SNPs.

**Table 6.** Concordance between ADMIXTURE and DAPC cluster assignments. Cross-tabulation of dominant ancestry assignments at K = 2 and K = 5.

**Table 7.** DAPC cluster sizes and eigenvalues. (A) Number of individuals assigned to each of the five clusters identified by K-means clustering (BIC criterion). (B) Eigenvalues and proportion of between-group variance explained by the four discriminant axes.

**Table 8.** Summary statistics of the centered IBS kinship matrix (PLINK `--make-rel`) computed from 26,474 QC-passing SNPs across 379 rice accessions. Thresholds of 0.125 and 0.25 correspond to second-degree and first-degree relatives, respectively.

**Figure 4.** ADMIXTURE analysis of 379 rice accessions. (a) Cross-validation error as a function of K (K = 1–12). The red point highlights K = 5. (b) Ancestry proportions at K = 2. (c) Ancestry proportions at K = 5. Each vertical bar represents one individual partitioned into colored segments proportional to estimated ancestry coefficients.

**Figure 5.** Discriminant Analysis of Principal Components (DAPC) of 379 rice accessions. (a) Scatter plot at K = 2. (b) Scatter plot at K = 5 showing five genetic clusters. (c) Posterior membership probabilities at K = 5.

**Figure 6.** Kinship matrix heatmap of 379 rice accessions based on 26,474 QC-passing SNPs. The centered IBS kinship matrix is visualized with hierarchical clustering (Euclidean distance, complete linkage).

**Supplementary Figure S1.** ADMIXTURE ancestry proportions at K = 12 for 379 rice accessions.
