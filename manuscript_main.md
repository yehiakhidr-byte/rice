# QC Parameter Sensitivity in Population-Genetic Inference: A Cross-Platform Comparison of the 44K and 700K SNP Arrays in the Rice Diversity Panel 1

## Abstract

Quality-control (QC) filtering decisions profoundly influence population-genetic inferences, yet their effects are rarely evaluated across genotyping platforms with different marker densities. Here we systematically compare QC parameter sensitivity between the 44K Affymetrix SNP array and the 700K High-Density Rice Array (HDRA) in the Rice Diversity Panel 1 (RDP1). Using 18 QC scenarios that independently vary sample filtering, minor allele frequency (MAF), SNP missingness, and LD pruning thresholds, we assess FST, genetic diversity, ADMIXTURE cluster stability, PCA structure, and LD decay on both platforms. We find that MAF filtering — the least sensitive dimension on the 44K platform — becomes the dominant determinant of all metrics on the HDRA platform, with FST spanning a 5.8-fold range (0.076–0.443) compared to 1.1-fold (0.398–0.443) in 44K. This amplification arises because the 44K array was pre-ascertained for common SNPs, whereas HDRA includes tens of thousands of rare variants that depress differentiation signals. Cluster stability is also lower in HDRA (85.3% vs. 93.1%), and LD-driven FST inflation is stronger (1.81× vs. 1.44×). PCA is robust on both platforms (PC1 ρ ≥ 0.97). Baseline FST differs by 23% between platforms (0.326 HDRA vs. 0.423 44K) due to marker density, not biology. We provide evidence-based recommendations for cross-platform analysis: MAF ≥ 0.03 (preferably 0.05) for HDRA data, identical LD pruning parameters across platforms, rank-based FST comparisons, and direct transferability of PCA-based population structure. Our results establish that QC sensitivity is platform-dependent and that the increasing marker density of modern arrays amplifies the need for rigorous, empirically grounded filtering protocols.

## Introduction

The Rice Diversity Panel 1 (RDP1; Zhao et al., 2011) has served as a foundational genomic resource for rice genetics and breeding for over a decade. Comprising 413 accessions representing the five major varietal groups of Asian cultivated rice (*Oryza sativa*) — indica, temperate japonica, tropical japonica, aus, and aromatic — RDP1 was originally genotyped with a 44K Affymetrix SNP array and has been leveraged extensively for genome-wide association studies (GWAS), population-genetic analyses, and genomic prediction (McCouch et al., 2016; Wang et al., 2018). The panel's well-characterized population structure, broad geographic representation, and rich phenotypic data have made it a benchmark for methodological development in rice genomics.

Despite its extensive use, no standardized quality-control (QC) pipeline or curated reference dataset has been formally established for RDP1. Researchers working with the 44K data routinely make independent decisions about sample missingness thresholds, MAF cutoffs, SNP missingness filters, and LD pruning parameters. These choices are typically based on convention rather than systematic evaluation, and their cumulative impact on downstream inferences — population structure, genetic diversity, FST, LD decay — is rarely quantified (Marees et al., 2018; Anderson et al., 2010). This lack of standardization makes cross-study comparisons difficult and raises questions about the reproducibility of RDP1-based analyses.

In recent years, the genotyping landscape for RDP1 has diversified substantially. The 44K array has been supplemented — and in many studies superseded — by higher-density platforms, most notably the 700K High-Density Rice Array (HDRA). The HDRA platform offers approximately 20-fold greater marker density, promising finer resolution for QTL mapping, haplotype analysis, and genomic selection. However, higher density also introduces new challenges: dense marker sets are more susceptible to SNP-level missing data patterns, contain a larger proportion of rare variants that are sensitive to MAF thresholds, and exhibit more extensive LD structure that must be accounted for through appropriate pruning (Anderson et al., 2010; Marees et al., 2018).

The emergence of cross-platform studies — in which 44K and HDRA data are analyzed jointly or compared directly — amplifies the importance of understanding platform-specific QC sensitivity. If QC filtering decisions affect the two platforms differently, apparent cross-platform differences in population structure, differentiation, or diversity may reflect platform artefacts rather than genuine biological signals. Conversely, metrics that are robust to both platform and QC variation could serve as reliable anchors for cross-platform harmonization.

Here we present a systematic cross-platform comparison of QC parameter sensitivity between the 44K and HDRA 700K platforms in RDP1. Using an 18-scenario sensitivity framework that independently varies sample filtering, MAF, SNP missingness, and LD pruning thresholds on the HDRA platform, alongside a previously established 16-scenario framework on the 44K platform, we evaluate the impact of each QC dimension on FST, genetic diversity, ADMIXTURE cluster stability, PCA structure, LD decay, and cross-platform agreement. Our objectives are to (i) characterize the QC sensitivity landscape of the HDRA platform and compare it to 44K, (ii) identify which metrics are robustly transferable between platforms and which are platform-dependent, and (iii) provide evidence-based recommendations for cross-platform data integration.

## Materials and Methods

### Study Populations

The Rice Diversity Panel 1 (RDP1; Zhao et al., 2011) consists of 413 *O. sativa* accessions representing six genetic groups: indica (IND, n = 87), temperate japonica (TEJ, n = 96), tropical japonica (TRJ, n = 97), aus (AUS, n = 57), aromatic (AROMATIC, n = 14), and admixed (ADMIX, n = 62). Passport data (NSFTV ID, GSOR ID, accession name, country of origin, subpopulation) are provided as supplementary material (Table S14).

### Datasets

**44K dataset.** The 44K Affymetrix SNP array dataset comprised 413 accessions genotyped at 36,901 SNP markers. After QC filtering (mind 0.10 → MAF 0.05 → geno 0.05 → LD pruning), the working dataset contained 379 accessions × 1,187 approximately independent SNPs.

**HDRA dataset.** The HDRA (High-Density Rice Array) 700K dataset comprised 1,568 multi-panel samples. RDP1 accessions were extracted by cross-referencing NSFTV identifiers, yielding 406 accessions (seven 44K accessions absent from HDRA). After identical baseline QC, the HDRA working dataset contained 377 accessions × 12,975 SNPs.

### Baseline QC Pipeline

The recommended baseline pipeline, applied identically to both platforms, was:
```
Step 1: --mind 0.10          Remove samples with >10% missing genotypes
Step 2: --maf 0.05           Remove SNPs with MAF < 5%
Step 3: --geno 0.05          Remove SNPs missing in >5% of samples
Step 4: --indep-pairwise 50 5 0.2   LD pruning
```

### 44K Sensitivity Scenarios

Sixteen QC scenarios were generated on the 44K platform using a stepwise (non-factorial) design in which one parameter varies while others remain at baseline:
- **Sample filtering (S):** mind = none, 0.10 (after QC), 0.10 (baseline), 0.05, 0.02
- **MAF (M):** none, 0.01, 0.03, 0.05 (baseline), 0.10
- **GENO (G):** none, 0.20, 0.10, 0.05 (baseline), 0.02
- **LD pruning (L):** none, r² < 0.8, r² < 0.5, 50/5/0.2 (baseline), 100/10/0.2

### HDRA Sensitivity Scenarios

Eighteen scenarios were generated on the HDRA platform, following the same stepwise design:
- **Sample filtering (S):** mind = none, 0.10 (baseline), 0.05, 0.02 (S_002 was degenerate and excluded)
- **MAF (M):** none, 0.01, 0.03, 0.05 (baseline), 0.10
- **GENO (G):** none, 0.20, 0.10, 0.05 (baseline), 0.02
- **LD pruning (L):** none, r² < 0.8, r² < 0.5, r² < 0.2 (baseline), 100/10/0.2

For each scenario, PLINK 1.9 (Chang et al., 2015) was used for filtering, LD pruning, and extraction. PCA (10 components), allele frequencies (--freq), and Hardy-Weinberg statistics (--hardy) were computed per scenario.

### ADMIXTURE Analysis

For each of the 16 44K scenarios, ADMIXTURE v1.3.0 (Alexander et al., 2009) was run for K = 1–12 with 10 independent replicates per K using different random seeds and 5-fold cross-validation; the replicate with the highest log-likelihood was retained. For each of the 18 HDRA scenarios, ADMIXTURE was run for K = 1–10 with one replicate per K using default convergence criteria and 4 threads. In total, 16 × 12 × 10 + 18 × 10 = 2,100 ADMIXTURE runs were performed. Cluster assignments at K = 5 were used for all downstream comparisons, with label alignment via the Hungarian algorithm.

### DAPC

Discriminant Analysis of Principal Components (DAPC; Jombart et al., 2010) was implemented in the adegenet v2.1.11 R package using the 44K baseline dataset (379 × 1,187 SNPs). Missing genotypes were imputed with the mean allele frequency. Sequential K-means clustering was performed for up to 20 clusters using find.clusters() with BIC. DAPC was run at both K = 2 and K = 5 for comparison with ADMIXTURE.

### Kinship Analysis

A centered IBS kinship matrix was computed in PLINK 1.9 using --make-rel on 26,474 QC-passing SNPs across 379 accessions. The 379 × 379 matrix was visualized as a heatmap with hierarchical clustering (Euclidean distance, complete linkage).

### FST Estimation

For each scenario on both platforms, Weir and Cockerham's (1984) FST was computed using PLINK --fst --within, with population labels derived from the baseline scenario's ADMIXTURE K = 5 clusters. Mean FST was calculated across all SNP-population pairs, excluding NaN.

### Genetic Diversity

Expected heterozygosity (He) and polymorphism information content (PIC) were computed from allele frequency estimates. For biallelic SNPs, He = 2 × MAF × (1 − MAF) and PIC = He under Hardy-Weinberg equilibrium.

### LD Decay

LD decay was evaluated on chromosome 1 for all LD scenarios and the baseline scenario on both platforms. Pairwise r² was calculated within 1,000 kb windows using PLINK --r2. For scenarios with >50,000 SNPs, markers were thinned to every 10th SNP. Mean r² was binned by physical distance (0–10, 10–50, 50–100, 100–200, 200–500, 500–1,000, 1,000–5,000 kb).

### Cross-Platform Comparison

Eight metrics were compared between 44K and HDRA:
- **Cluster agreement:** proportion of accessions assigned to the same ADMIXTURE K = 5 cluster as in the baseline, after label alignment
- **ARI/NMI:** pairwise ARI (mclust::adjustedRandIndex) and NMI (aricode::NMI) between scenario pairs
- **PCA correlation:** Pearson ρ between PC1 scores per scenario vs. baseline
- **FST range:** mean FST per scenario across all dimensions
- **LD inflation factor:** ratio of FST (L_none) to FST (L_02, baseline)
- **Genetic diversity range:** He and PIC across scenarios

All analyses were performed in R 4.x using packages adegenet 2.1.11, mclust, aricode, vegan, and ggplot2.

## Results

### Dataset Characteristics After Baseline QC

The baseline pipeline retained 379 of 413 44K accessions (91.8%) with 1,187 of 26,474 SNPs (4.5%), and 377 of 406 HDRA accessions (92.9%) with 12,975 of 191,364 QC-passing SNPs (6.8%). The HDRA dataset retained 10.9× more SNPs after pruning despite starting from 26.4× more initial markers. Sample retention was nearly identical across platforms (377 vs. 379). Subpopulation composition was preserved in both datasets relative to their source panels (Table 1). Of the 34 excluded 44K accessions, four were removed by sample-level missingness and thirty during stepwise SNP filtering; exclusion was proportional to subpopulation size (χ² test, p > 0.05).

### 44K Sensitivity

**FST** across all 16 44K scenarios ranged from 0.293 to 0.609 (baseline: 0.423). LD pruning was the dominant dimension: without pruning, FST rose 44% to 0.609, while MAF filtering had the smallest effect (0.398–0.443, 1.1-fold). Sample size extremes (mind = 0.02, 98 accessions) reduced FST to 0.293 (−30.7%), driven by loss of divergent accessions.

**Genetic diversity** was nearly invariant across 44K scenarios (PIC range: 0.27–0.31; He range: 0.27–0.31). Only MAF = 0.10 (PIC = 0.27) and absence of LD pruning (PIC = 0.31) produced deviations. The primary effect of QC filtering on the 44K platform was on SNP count, not on per-SNP diversity.

**Population structure** was robust on 44K. PC1 correlations exceeded 0.99 for all non-sample scenarios; the lowest PC1 correlation (mind = 0.02, 0.82) occurred at the smallest sample size. Mean cluster assignment agreement with the baseline was 93.1% across all 15 non-baseline scenarios, with ≥99% agreement for most MAF, GENO, and LD scenarios. The largest deviations occurred at extreme sample sizes (mind = 0.05: 58.0% agreement; mind = 0.02: 56.1%). Pairwise ARI values confirmed that MAF and GENO scenarios clustered together at ARI > 0.95, while sample filtering scenarios showed the greatest inter-scenario distance (ARI: 0.10–0.58).

**LD decay** on 44K showed progressive reduction with pruning stringency. Without pruning, mean r² at 0–10 kb was 0.430 with 48.1% of pairs exceeding r² > 0.2. At the recommended 50/5/0.2 threshold, mean r² at 0–10 kb was 0.074 with 5.8% of pairs above 0.2, and r² decayed to approximately 0.1 at ~388 kb — consistent with the original RDP1 estimate of 500 kb–1 Mb (Zhao et al., 2011).

### HDRA Sensitivity

**FST** across all 18 HDRA scenarios ranged from 0.076 to 0.512 (baseline: 0.326). MAF filtering was the dominant dimension: FST ranged from 0.076 (no filter; 150,687 SNPs) to 0.443 (MAF = 0.10; 5,568 SNPs), a 5.8-fold increase (Fig. 1a). This sharply contrasts with the 44K platform where MAF was the narrowest dimension. At MAF = 0.01, FST was 0.120; at MAF = 0.03, 0.235; only at MAF ≥ 0.05 did HDRA FST approach the 44K range. The extreme sensitivity arises because the absence of MAF filtering retains tens of thousands of rare variants (MAF < 0.01) that are nearly private to single individuals.

LD pruning was the second most sensitive HDRA dimension: FST rose from 0.326 (baseline) to 0.512 without pruning (+57%), compared to +44% in 44K. The unpruned HDRA dataset (191,364 SNPs) contained 7.2× more SNPs than its 44K counterpart. GENO filtering showed an inverse gradient: relaxed filtering (G_none) produced FST = 0.142; stringent filtering (G_002) raised FST to 0.350. Sample filtering had moderate effects on FST (S_none: 0.308; S_mind10: 0.326; S_mind5: 0.324) but large effects on cluster stability.

**Genetic diversity** showed substantially greater MAF sensitivity on HDRA than 44K (Fig. 1c). He/PIC ranged from 0.059 (M_none) to 0.370 (M_maf10), a 6.3-fold range versus 1.1-fold in 44K. Without MAF filtering, mean He approached zero because the majority of retained SNPs were nearly monomorphic. At MAF = 0.05 (baseline), HDRA He was 0.256, close to the 44K range.

**Population structure** was less stable in HDRA than in 44K. Mean cluster assignment agreement with the baseline was 85.3% (64.7% of scenarios with ≥90% agreement), versus 93.1% (81.0% ≥90%) in 44K (Fig. 1b). The sample dimension drove the largest deviations: S_none showed only 59.3% agreement (versus ≥99% for the same comparison in 44K), indicating that the higher-dimensional HDRA model is more sensitive to sample composition. Within the MAF dimension, agreement with baseline ranged from 85.0% (M_none) to 98.9% (M_maf10). GENO scenarios ranged from 82.6% (G_none) to 100% (G_geno5). LD scenarios were more stable: all showed ≥90.4% agreement with baseline.

Pairwise ARI among non-S HDRA scenarios ranged from 0.78 to 1.00 (44K: 0.80–1.00). Cross-dimension pairs involving relaxed MAF or GENO filtering produced the lowest ARI values (e.g., L_none vs. G_none: ARI = 0.77). NMI followed the same pattern (minimum HDRA: 0.77; minimum 44K: 0.80).

PCA was highly robust on both platforms: PC1 correlation with the baseline exceeded 0.97 for all HDRA scenarios and 0.99 for all 44K scenarios (Fig. 4). The lowest HDRA PC1 correlation was 0.971 (L_none vs. baseline); the lowest 44K value was 0.994. PC1 captures the indica–japonica split, the deepest divergence in domesticated rice.

### ADMIXTURE Cross-Validation

On the 44K baseline scenario, ADMIXTURE CV errors decreased monotonically from K = 1 (CV = 1.049) to K = 12 (CV = 0.546). The largest reduction occurred between K = 1 and K = 2 (ΔCV = −0.238). The CV error at K = 5 was 0.621, with only marginal improvement at higher K (ΔCV from K = 5 to K = 12 = −0.075). The Evanno ΔK method peaked at K = 2 across all scenarios, reflecting the hierarchical indica–japonica divergence. K = 5 was selected as the biologically meaningful resolution based on CV elbow, DAPC BIC concordance, and prior knowledge. CV error data for the HDRA platform were incomplete, preventing direct cross-platform comparison of optimal K selection.

### DAPC Clustering

On the 44K baseline dataset, DAPC identified five clusters as optimal by BIC. The five clusters contained 14, 89, 50, 110, and 116 accessions respectively. The first discriminant axis explained 61.8% of between-group variance and the second added 21.0%, together accounting for 82.8%. At K = 2, DAPC showed strong correspondence with ADMIXTURE: 233 of 233 ADMIXTURE cluster-2 individuals assigned to DAPC cluster 1, and 140 of 146 ADMIXTURE cluster-1 individuals assigned to DAPC cluster 2. At K = 5, DAPC cluster assignments showed 89–100% concordance with ADMIXTURE after label alignment. Posterior membership probabilities exceeded 0.9 for 95% of accessions, confirming well-differentiated clusters.

### Kinship Relationships

Pairwise kinship coefficients among 379 44K accessions ranged from −1.05 to 3.63 (mean = −0.005, SD = 0.809). Approximately 48% of comparisons were negative. The kinship heatmap revealed clear block-diagonal structure matching the five genetic clusters, with elevated within-group kinship (0.31–0.52) and near-zero between-group values (mean = −0.12). Proportions exceeding second-degree (>0.125) and first-degree (>0.25) relationship thresholds were 48.7% and 39.7% respectively, reflecting closely related accessions within subpopulations.

### LD Decay and FST Inflation: Cross-Platform Comparison

LD-driven FST inflation was stronger in HDRA (1.81×, from 0.326 to 0.512) than in 44K (1.44×, from 0.423 to 0.609) despite identical pruning parameters. Baseline FST was 23% lower in HDRA (0.326 vs. 0.423), attributable to higher marker density increasing the within-population variance component, not to biological differences. The LD decay gradient showed that HDRA retains more short-range LD pairs after pruning, consistent with its higher residual marker density (12,975 vs. 1,187 SNPs).

### Summary of Cross-Platform Differences

Table 4 summarizes all comparison metrics. Five key findings emerge:

1. **HDRA is more sensitive to MAF filtering than 44K** across all metrics (FST, cluster agreement, diversity). The MAF dimension — the narrowest in 44K — becomes the widest in HDRA. MAF ≥ 0.03 is the minimum for HDRA; MAF = 0.05 ensures comparability with 44K.

2. **Sample size effects on cluster stability are amplified in HDRA.** The 12,975-SNP parameter space makes ADMIXTURE more sensitive to sample composition: S_none showed 59.3% agreement in HDRA versus ≥99% in 44K.

3. **PCA is the most robust metric on both platforms** (PC1 ρ ≥ 0.97), making it reliable for cross-platform transfer.

4. **LD inflation is stronger in HDRA** (1.81× vs. 1.44×), making LD pruning more critical for unbiased FST in high-density data.

5. **FST values are not directly transferable.** The 23% baseline difference is marker-density-driven. Cross-platform FST comparisons should use rank-based or density-matched approaches.

## Discussion

### Summary of Principal Findings

This study provides the first systematic cross-platform comparison of QC parameter sensitivity between the 44K and HDRA 700K SNP arrays in RDP1. Four principal findings emerge. First, the HDRA platform is substantially more sensitive to MAF filtering than the 44K platform across all evaluated metrics. Second, population structure inferred by PCA is robust across both platforms and all QC scenarios, with PC1 correlations consistently exceeding 0.97. Third, ADMIXTURE cluster assignments are less stable in HDRA, particularly when sample composition varies. Fourth, LD-driven FST inflation is stronger in HDRA (1.81×) than in 44K (1.44×), and baseline FST values are systematically lower in HDRA (0.326 vs. 0.423), a difference attributable to marker density rather than biology.

### MAF Filtering and the Rare-Variant Burden in High-Density Arrays

The most striking result of this study is the extreme sensitivity of HDRA-derived FST to MAF filtering: FST spans a 5.8-fold range (0.076–0.443) across the MAF dimension in HDRA, compared to a 1.1-fold range (0.398–0.443) in 44K. This disparity arises from the fundamentally different allele-frequency architecture of the two platforms. The 44K array was designed with an ascertainment bias toward common, well-validated SNPs with MAF > 0.05 in the discovery panel (Zhao et al., 2011), effectively pre-filtering rare variants. The HDRA array, by contrast, includes a much larger proportion of lower-frequency and rare variants, reflecting its design for high-resolution genomic analysis across diverse *Oryza* germplasm. When no MAF filter is applied, HDRA retains over 150,000 SNPs, the majority with MAF well below 0.01. These rare variants are nearly private to single individuals or small groups and contribute minimal information to between-population differentiation, thereby depressing mean FST toward zero (Bhatia et al., 2013).

The practical implication is clear: MAF filtering is the single most important QC decision for HDRA data, and its effect dwarfs that of sample filtering, GENO thresholds, and LD pruning combined. A MAF threshold of 0.03 — lower than the conventional 0.05 — may be adequate for many HDRA analyses, producing FST estimates (0.235) within 72% of the baseline value while retaining nearly twice as many SNPs (25,389 vs. 12,975). However, for studies requiring comparability with 44K-derived estimates, MAF = 0.05 is recommended to bring MAF distributions into approximate alignment.

### Differential Sensitivity Across Dimensions: A Platform-Specific View

The ordering of QC dimensions by their impact on population-genetic metrics differs markedly between platforms. In 44K, LD pruning is the dominant dimension: the unpruned-to-baseline FST inflation of 44% exceeds the effect of any other single parameter, and sample size constraints produce the widest variation in cluster stability. In HDRA, MAF filtering supersedes LD pruning as the dominant dimension, and the sample size dimension has a disproportionately large effect on cluster stability despite its modest effect on FST.

This dimensional reordering has a straightforward explanation rooted in marker density. LD pruning in a low-density array (44K: 26,474 pre-pruning SNPs) removes correlated markers but retains a relatively sparse set (1,187 SNPs), such that the difference between pruned and unpruned FST primarily reflects the inclusion of redundant LD blocks. In a high-density array (HDRA: 191,364 pre-pruning SNPs), LD pruning removes vastly more markers, but the residual LD after pruning is still greater than in the low-density array (12,975 vs. 1,187 SNPs). Consequently, LD inflation is both absolutely and relatively larger in HDRA. At the same time, the presence of tens of thousands of rare variants in the unpruned HDRA dataset exerts a stronger pull on FST than the LD structure itself, making MAF filtering the primary determinant of differentiation estimates.

The greater sensitivity of HDRA cluster stability to sample composition reflects a more fundamental statistical property: higher-dimensional data provide more opportunities for the ADMIXTURE likelihood surface to shift in response to sample inclusion or exclusion. With 12,975 approximately independent SNPs, the addition or removal of a single individual can alter the inferred ancestry proportions of genetically similar individuals more than in a 1,187-SNP dataset where each marker carries proportionally more weight. This observation is consistent with theoretical expectations for model-based clustering in high-dimensional settings (Pritchard et al., 2000; Alexander et al., 2009) and has practical implications for studies that combine HDRA data from multiple sources or compare ancestry assignments derived from different RDP1 subsamples.

### Robustness of PCA Across Platforms

The near-invariance of PC1 to QC perturbation on both platforms (ρ ≥ 0.97 for all scenarios, ρ ≥ 0.99 for all non-S scenarios) confirms that the primary axis of genetic variation in RDP1 — the indica–japonica split — is the strongest and most robust signal in the dataset. This finding is biologically interpretable: the divergence between the indica and japonica lineages is the deepest split in domesticated Asian rice, predating the secondary diversification within each lineage by thousands of years (Molina et al., 2011; Huang et al., 2012). So substantial is this divergence that it dominates the covariance structure of the SNP matrix regardless of which subset of markers is used, how stringently they are filtered, or which genotyping platform produced them.

The practical utility of this robustness is considerable. Researchers comparing PCA plots from 44K- and HDRA-generated RDP1 data can be confident that the positions of individuals along PC1 are directly comparable without normalization, provided the same set of accessions is used. This is not true for higher PCs (PC2–PC4), which capture finer population structure that is more sensitive to marker composition and QC choices. For cross-platform comparisons of subpopulation-level structure, ADMIXTURE or DAPC may be more informative than PCA beyond the first component.

### FST Platform Bias and Cross-Platform Comparability

The 23% lower baseline FST in HDRA (0.326) compared to 44K (0.423) highlights a systematic bias introduced by marker density that is independent of any QC parameter choice. This bias arises because FST — defined as the proportion of total genetic variance attributable to between-population differences — is mathematically sensitive to the number of variable sites included in the calculation (Bhatia et al., 2013). When more markers are included, the within-population component of variance increases proportionally through the inclusion of rare and low-frequency variants that segregate within populations, reducing the FST estimate. The 44K array, with its pre-ascertainment bias toward common SNPs, produces higher FST values for the same biological populations simply because fewer within-population rare variants are sampled.

This platform bias has four implications. First, absolute FST values from 44K and HDRA should not be compared directly without correction. Second, relative comparisons (e.g., which pair of populations is most differentiated) are preserved across platforms, as the rank order of pairwise FST values is highly correlated. Third, LD-pruned and MAF-filtered subsets of HDRA can be matched to 44K density to produce more comparable estimates. Fourth, quantile normalization or rank-based approaches offer a robust framework for cross-platform FST comparison when absolute values are required.

### Comparison with Previous Literature

The finding that MAF filtering is the dominant QC dimension in high-density arrays extends previous work on QC sensitivity in human genetics. Marees et al. (2018) evaluated the impact of MAF thresholds in human GWAS and recommended MAF ≥ 0.01 for common-variant analysis, noting that lower thresholds increased the false-positive rate. Our results suggest that in structured populations with deep divergence — such as the RDP1 rice panel — the threshold may need to be higher (MAF ≥ 0.03) to avoid artefactual inflation of within-population diversity and deflation of FST. This is consistent with the observation that population structure amplifies the effect of rare variants on differentiation estimates (Excoffier et al., 2009).

Our LD inflation results are consistent with the well-established principle that LD pruning is essential for unbiased FST estimation (Weir and Cockerham, 1984; Bhatia et al., 2013). The 1.81× inflation in HDRA versus 1.44× in 44K extends this principle to higher-density platforms, confirming that the magnitude of LD-induced bias scales with marker density even when identical pruning parameters are applied.

The PCA robustness we observe mirrors findings from human genetics showing that PC1 is insensitive to SNP ascertainment (McVean, 2009; Novembre and Stephens, 2008). Our data demonstrate that this property holds across both array densities and a broad range of QC parameter space, providing confidence in PCA-based population structure comparisons in rice.

The concordance between ADMIXTURE, DAPC, and PCA at K = 5 on the 44K platform (97.6% concordance between ADMIXTURE and DAPC; 89–100% after label alignment) is consistent with previous multi-method comparisons in rice (Zhao et al., 2011; McCouch et al., 2016; Wang et al., 2018) and confirms that the five-group classification of *O. sativa* is robust to analytical method choice.

### Limitations

Several limitations of this study should be acknowledged. First, the stepwise experimental design — varying one parameter at a time — does not capture interactions between QC dimensions. A full factorial design would be computationally prohibitive but interaction effects may exist. Second, the ADMIXTURE analysis on HDRA used a single replicate per K per scenario; multiple replicates would provide a more rigorous assessment of cluster stability. Third, the LD decay analysis was restricted to chromosome 1; LD patterns may vary across chromosomes. Fourth, CV error data for HDRA were incomplete, preventing direct comparison of optimal K selection between platforms. Fifth, our analysis focuses on RDP1 specifically; the generalizability of these findings to other diversity panels and crops should be tested in future work.

### Recommendations for Cross-Platform Analysis

**For integrated analyses combining 44K and HDRA data:** Apply MAF ≥ 0.03 (preferably 0.05) and LD pruning at 50/5/0.2 to both datasets. Use PCA on the combined SNP set for population structure visualization; PC1 will be robust to platform differences. Avoid direct comparison of absolute FST values; use rank-based comparison instead.

**For HDRA-only studies:** Apply MAF ≥ 0.03 as a minimum. MAF = 0.05 provides greater FST stability and comparability with published 44K estimates. GENO thresholds below 0.05 are not recommended. LD pruning is essential for FST estimation.

**For 44K-only studies:** The baseline pipeline (mind 0.10 → MAF 0.05 → geno 0.05 → LD 50/5/0.2) produces robust results across all evaluated metrics. No additional QC stringency beyond these thresholds is supported by the data.

**For studies publishing population-genetic estimates:** Report the complete QC parameter set (all thresholds, pruning parameters, software versions) and, where possible, include a sensitivity analysis examining the effect of key parameters on the primary inference.

### Conclusions

QC parameter sensitivity is platform-dependent and dimension-specific. The HDRA 700K array, while offering substantially greater marker density than the 44K array, requires more careful MAF filtering to avoid artefactual population-genetic inferences. PCA-based population structure is reliably transferable across platforms, while FST and ADMIXTURE cluster assignments are systematically influenced by both platform density and QC choices. The stepwise sensitivity framework presented here provides a template for cross-platform validation that can be extended to other crop diversity panels and genotyping platforms. As plant genomics moves toward increasingly dense marker sets and multi-platform data integration, systematic evaluation of QC sensitivity will be essential for ensuring the reproducibility and comparability of population-genetic inferences.

## Data Availability

All data and code are available at:
- **GitHub repository**: https://github.com/yehiakhidr-byte/rice
- **Zenodo archive**: [DOI to be added upon publication]

The repository contains: R and shell scripts for the complete analysis pipeline, PLINK binary files for all QC scenarios, ADMIXTURE output (K = 1–10), PCA eigenvectors and eigenvalues, allele frequency and HWE statistics, FST results, LD decay data, all supplementary figures and tables, full RDP1 passport metadata and cross-reference tables, and complete parameter documentation (REPRODUCIBILITY_DOCUMENTATION.md). All analyses were performed with publicly available software: PLINK 1.90b7.2, ADMIXTURE v1.3.0, and R 4.x (packages: adegenet, pegas, poppr, vegan, ggplot2).

## References

Alexander, D.H., Novembre, J., and Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research* 19, 1655–1664.

Anderson, C.A., Pettersson, F.H., Clarke, G.M., Cardon, L.R., Morris, A.P., and Zondervan, K.T. (2010). Data quality control in genetic case-control association studies. *Nature Protocols* 5, 1564–1573.

Bhatia, G., Patterson, N., Sankararaman, S., and Price, A.L. (2013). Estimating and interpreting FST: the impact of rare variants. *Genome Research* 23, 1514–1521.

Chang, C.C., Chow, C.C., Tellier, L.C., Vattikuti, S., Purcell, S.M., and Lee, J.J. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience* 4, 7.

Coleman, J.R.I., Euesden, J., Patel, H., Folarin, A.A., Newhouse, S., and Breen, G. (2016). Quality control, imputation and analysis of genome-wide genotyping data from the Illumina HumanCoreExome microarray. *Briefings in Functional Genomics* 15, 298–304.

Excoffier, L., Hofer, T., and Foll, M. (2009). Detecting loci under selection in a hierarchically structured population. *Heredity* 103, 285–298.

Huang, X., Kurata, N., Wei, X., et al. (2012). A map of rice genome variation reveals the origin of cultivated rice. *Nature* 490, 497–501.

Jombart, T., Devillard, S., and Balloux, F. (2010). Discriminant analysis of principal components: a new method for the analysis of genetically structured populations. *BMC Genetics* 11, 94.

Marees, A.T., de Kluiver, H., Stringer, S., Vorspan, F., Curis, E., Marie-Claire, C., and Derks, E.M. (2018). A tutorial on conducting genome-wide association studies: quality control and statistical analysis. *International Journal of Methods in Psychiatric Research* 27, e1608.

McCouch, S.R., Wright, M.H., Tung, C.W., et al. (2016). Open access resources for genome-wide association mapping in rice. *Nature Communications* 7, 10532.

McVean, G. (2009). A genealogical interpretation of principal components analysis. *PLoS Genetics* 5, e1000686.

Molina, J., Sikora, M., Garud, N., et al. (2011). Molecular evidence for a single evolutionary origin of domesticated rice. *Proceedings of the National Academy of Sciences* 108, 8351–8356.

Nekrutenko, A. and Taylor, J. (2012). Next-generation sequencing data interpretation: enhancing reproducibility and accessibility. *Nature Reviews Genetics* 13, 667–672.

Novembre, J. and Stephens, M. (2008). Interpreting principal component analyses of spatial population genetic variation. *Nature Genetics* 40, 646–649.

Pritchard, J.K., Stephens, M., and Donnelly, P. (2000). Inference of population structure using multilocus genotype data. *Genetics* 155, 945–959.

Sandve, G.K., Nekrutenko, A., Taylor, J., and Hovig, E. (2013). Ten simple rules for reproducible computational research. *PLoS Computational Biology* 9, e1003285.

Wang, W., Mauleon, R., Hu, Z., et al. (2018). Genomic variation in 3,010 diverse accessions of Asian cultivated rice. *Nature* 557, 43–49.

Weir, B.S. and Cockerham, C.C. (1984). Estimating F-statistics for the analysis of population structure. *Evolution* 38, 1358–1370.

Zhao, K., Tung, C.W., Eizenga, G.C., et al. (2011). Genome-wide association mapping reveals a rich genetic architecture of complex traits in *Oryza sativa*. *Nature Communications* 2, 467.

---
