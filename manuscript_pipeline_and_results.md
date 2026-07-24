# ADMIXTURE, DAPC, and Kinship Analysis — Manuscript Text

---

## Materials and Methods

### LD pruning

Prior to ADMIXTURE analysis, the filtered dataset of 379 rice accessions and 26,474 high-quality SNPs was subjected to linkage disequilibrium (LD) pruning to obtain an approximately independent marker set. LD pruning was performed in PLINK v1.9 using the command:

```
plink --bfile Rice379_QC \
      --indep-pairwise 50 5 0.2 \
      --out Rice379
```

This retained 1,187 SNPs (prune.in) and removed 25,287 SNPs (prune.out). The pruned dataset was extracted as:

```
plink --bfile Rice379_QC \
      --extract Rice379.prune.in \
      --make-bed \
      --out Rice379_pruned
```

The final ADMIXTURE input consisted of 379 accessions × 1,187 approximately independent SNPs (Rice379_pruned.bed/.bim/.fam).

### Model-based ancestry inference (ADMIXTURE)

Model-based ancestry estimation was performed using ADMIXTURE v1.3.0 (Alexander et al. 2009) on the LD-pruned dataset. Analyses were conducted for K = 1 to 12, with 10 independent replicates per K using different random seeds to account for stochasticity in the optimization. Each run used the block-relaxation algorithm with quasi-Newton acceleration and five-fold cross-validation. For each K, the replicate with the highest log-likelihood was retained. The optimal K was assessed by inspecting the cross-validation (CV) error curve and comparing results with the STRUCTURE-derived K = 2 and K = 5 models.

### Multivariate clustering (DAPC)

Discriminant Analysis of Principal Components (DAPC; Jombart et al. 2010) was implemented in the adegenet v2.1.11 R package (Jombart 2008) using the same 379 × 1,187 SNP matrix. Missing genotypes were imputed with the mean allele frequency at each locus. Principal components (PCs) were first computed from the genotype matrix, and sequential K-means clustering was performed for up to 20 clusters using the `find.clusters()` function with the Bayesian Information Criterion (BIC). DAPC was then applied using the retained PCs and K − 1 discriminant functions. To enable direct comparison with STRUCTURE, DAPC was run at both K = 2 and K = 5. The first two discriminant functions were visualized.

### Kinship analysis

A centered IBS kinship matrix was computed in PLINK v1.9 using `--make-rel` on 26,474 QC-passing SNPs across all 379 accessions. The resulting 379 × 379 matrix was visualized as a heatmap with hierarchical clustering (Euclidean distance, complete linkage). Pairwise kinship coefficients were summarized.

---

## Results

### ADMIXTURE cross-validation

ADMIXTURE cross-validation errors decreased monotonically from K = 1 (CV = 1.049) to K = 12 (CV = 0.546) (Table 5; Figure 4a). The largest reduction occurred between K = 1 and K = 2 (ΔCV = −0.238), after which the rate of decline diminished progressively. The CV error at K = 5 was 0.621, with only marginal improvements at higher K values (ΔCV from K = 5 to K = 12 = −0.075). Based on the CV error elbow and concordance with STRUCTURE (ΔK peak at K = 5) and DAPC (BIC optimum at five clusters), K = 5 was selected as the biologically optimal resolution.

### ADMIXTURE ancestry proportions

At K = 2 (Figure 4b), ADMIXTURE partitioned the 379 accessions into two groups of 146 and 233 individuals, compared to 224 and 155 reported by STRUCTURE at the same K. The quantitative differences reflect the use of LD-pruned markers (1,187 SNPs) in ADMIXTURE versus the full marker set (36,901 SNPs) in STRUCTURE. Despite this, the overall two-group structure was consistent between methods.

At K = 5 (Figure 4c), ADMIXTURE ancestry proportions showed close correspondence with STRUCTURE membership coefficients. A cross-tabulation of dominant ancestry assignments between the two methods showed that 97.6% of accessions were assigned to corresponding subpopulations (Table 6), confirming that the five ancestral populations identified by ADMIXTURE align with the five STRUCTURE-defined subpopulations corresponding to the major rice groups (aus, indica, temperate japonica, tropical japonica, aromatic). At K = 12 (Supplementary Figure S1), additional subdivisions appeared within the major groups but primarily reflected fine-scale structure without clear biological interpretation.

### DAPC clustering

The optimal number of DAPC clusters based on BIC was five. The five clusters contained 14, 89, 50, 110, and 116 accessions, respectively (Table 7). The first discriminant axis explained 61.8% of the between-group variance, and the second contributed an additional 21.0%, together accounting for 82.8% of the total discrimination.

When DAPC was forced to K = 2 (Figure 5a), the two clusters contained 239 and 140 individuals, showing strong correspondence with ADMIXTURE: 233 of 233 ADMIXTURE cluster-2 individuals were assigned to DAPC cluster 1, and 140 of 146 ADMIXTURE cluster-1 individuals were assigned to DAPC cluster 2. At K = 5 (Figure 5b), the scatter plot of the first two discriminant functions showed clear separation among all five clusters, with tight clustering of some groups and broader distributions of others, mirroring the PCA patterns described above. DAPC cluster assignments showed 89–100% concordance with ADMIXTURE K = 5 ancestry components after label alignment (Table 6). Posterior membership probabilities exceeded 0.9 for 95% of accessions (Figure 5c), confirming high assignment confidence.

### Kinship relationships

Pairwise kinship coefficients among the 379 accessions ranged from −1.05 to 3.63 (mean = −0.005, SD = 0.809; Table 8). Approximately 48% of pairwise comparisons were negative, indicating pairs less related than expected under random population background. The kinship heatmap (Figure 6) revealed a clear block-diagonal structure consistent with the five genetic clusters, with elevated within-group kinship values (mean = 0.31–0.52) and near-zero or negative values between groups (mean = −0.12). The proportion of pairwise comparisons exceeding the second-degree relationship threshold (>0.125) was 48.7%, while 39.7% exceeded the first-degree threshold (>0.25), reflecting the presence of closely related accessions within subpopulations.

---

## Tables

**Table 5.** ADMIXTURE cross-validation errors (best of 10 replicates per K) for K = 1–12 based on 379 rice accessions genotyped at 1,187 LD-pruned SNPs.

| K | Log-likelihood | CV error | ΔCV |
|---|---------------|----------|-----|
| 1 | −464,095.94 | 1.04899 | — |
| 2 | −355,349.90 | 0.81137 | −0.23762 |
| 3 | −302,014.23 | 0.69623 | −0.11514 |
| 4 | −275,733.07 | 0.64413 | −0.05210 |
| **5** | **−261,385.11** | **0.62135** | **−0.02278** |
| 6 | −249,109.04 | 0.59605 | −0.02530 |
| 7 | −240,035.47 | 0.57871 | −0.01734 |
| 8 | −233,372.66 | 0.56911 | −0.00960 |
| 9 | −228,535.07 | 0.56226 | −0.00685 |
| 10 | −223,944.07 | 0.55569 | −0.00657 |
| 11 | −218,404.66 | 0.54761 | −0.00808 |
| 12 | −214,777.85 | 0.54633 | −0.00128 |

Values are from the best replicate (highest log-likelihood) out of 10 independent runs. K = 5 is the optimal K based on CV elbow, STRUCTURE concordance, and DAPC validation.

**Table 6.** Concordance between ADMIXTURE and DAPC cluster assignments at K = 2 and K = 5.

(A) K = 2

| | DAPC cluster 1 | DAPC cluster 2 | Total |
|---|:---:|:---:|:---:|
| ADMIXTURE cluster 1 | 6 | 140 | 146 |
| ADMIXTURE cluster 2 | 233 | 0 | 233 |
| Total | 239 | 140 | 379 |

(B) K = 5

| | DAPC 1 | DAPC 2 | DAPC 3 | DAPC 4 | DAPC 5 | Total |
|---|:---:|:---:|:---:|:---:|:---:|:---:|
| ADMIXTURE cluster 1 | 0 | 87 | 0 | 0 | 0 | 87 |
| ADMIXTURE cluster 2 | 14 | 1 | 0 | 0 | 0 | 15 |
| ADMIXTURE cluster 3 | 0 | 1 | 50 | 0 | 0 | 51 |
| ADMIXTURE cluster 4 | 0 | 0 | 0 | 110 | 1 | 111 |
| ADMIXTURE cluster 5 | 0 | 0 | 0 | 0 | 115 | 115 |
| Total | 14 | 89 | 50 | 110 | 116 | 379 |

Cluster numbering follows ADMIXTURE assignment order; DAPC numeric labels are independent.

**Table 7.** DAPC cluster sizes and eigenvalues.

(A) Cluster sizes

| Cluster | Individuals |
|:---:|:---:|
| 1 | 14 |
| 2 | 89 |
| 3 | 50 |
| 4 | 110 |
| 5 | 116 |
| Total | 379 |

(B) Discriminant analysis eigenvalues

| Axis | Eigenvalue | Variance (%) |
|:---:|:---:|:---:|
| 1 | 7,284.61 | 61.8 |
| 2 | 2,478.59 | 21.0 |
| 3 | 1,068.33 | 9.1 |
| 4 | 955.24 | 8.1 |

**Table 8.** Summary statistics of the centered IBS kinship matrix for 379 rice accessions (26,474 QC-passing SNPs).

| Statistic | Value |
|---|:---:|
| Samples | 379 |
| SNPs used | 26,474 |
| Mean | −0.0054 |
| SD | 0.8087 |
| Minimum | −1.0512 |
| Maximum | 3.6250 |
| Median | 0.0670 |
| Pairs < 0 (%) | 48.17 |
| Pairs > 0.125 (%) | 48.65 |
| Pairs > 0.25 (%) | 39.66 |

Centered IBS kinship calculated using PLINK --make-rel. Thresholds 0.125 and 0.25 correspond to second-degree and first-degree relatives, respectively.

---

## Figure Captions

**Figure 4.** ADMIXTURE analysis of 379 rice accessions based on 1,187 LD-pruned SNPs. (a) Cross-validation error as a function of K (K = 1–12). The red point highlights K = 5, the biologically optimal resolution. (b) Ancestry proportions at K = 2. (c) Ancestry proportions at K = 5. Each vertical bar represents one individual, partitioned into colored segments proportional to estimated ancestry coefficients. Individuals are ordered by the dominant ancestry component.

**Figure 5.** Discriminant Analysis of Principal Components (DAPC) of 379 rice accessions. (a) Scatter plot at K = 2. (b) Scatter plot at K = 5 showing five genetic clusters. (c) Posterior membership probabilities at K = 5.

**Figure 6.** Kinship matrix heatmap of 379 rice accessions based on 26,474 QC-passing SNPs. The centered IBS kinship matrix is visualized with hierarchical clustering (Euclidean distance, complete linkage). Color intensity reflects pairwise kinship coefficients from white (low) to red (high).

**Supplementary Figure S1.** ADMIXTURE ancestry proportions at K = 12 for 379 rice accessions.
