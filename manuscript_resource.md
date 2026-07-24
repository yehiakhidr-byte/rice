# RiceDP1-QC: A Curated Genomic Resource and Reproducible Quality-Control Workflow for the Rice Diversity Panel 1

## Background & Context

The Rice Diversity Panel 1 (RDP1; Zhao et al., 2011) is one of the most widely used genomic resources in rice genetics and breeding. Comprising 413 accessions of *Oryza sativa* genotyped with a 44K SNP array, it has been the foundation for genome-wide association studies (GWAS), population-genetic analyses, and genomic prediction in rice for over a decade. Despite its extensive use, no standardized quality-control (QC) pipeline or curated reference dataset has been formally established for RDP1.

Researchers working with RDP1 routinely make independent decisions about sample missingness thresholds, minor allele frequency (MAF) cutoffs, SNP missingness filters, and linkage disequilibrium (LD) pruning parameters. These choices are often based on convention rather than systematic evaluation, and their cumulative impact on downstream inferences—population structure, genetic diversity, FST, LD decay—is rarely quantified. This lack of standardization makes cross-study comparisons difficult and raises questions about the reproducibility of RDP1-based analyses.

Several recent reviews have emphasized the need for reproducible analytical workflows in population genomics (Nekrutenko & Taylor, 2012; Sandve et al., 2013), and the importance of documenting QC decisions has been recognized by journals and funding agencies. However, practical implementations for specific crop diversity panels remain scarce, and even widely used panels like RDP1 lack an accompanying validated pipeline.

Here we present RiceDP1-QC: a curated genomic resource for RDP1 that includes (1) a fully documented, stepwise QC pipeline with empirically validated threshold recommendations, (2) a set of 16 systematically generated QC scenarios with complete population-genomic analyses, (3) a comparison of population-structure inference methods (PCA, ADMIXTURE, DAPC) and their sensitivity to QC decisions, and (4) a reproducible computational workflow deposited in a public GitHub repository with a Zenodo DOI. The resource is designed to serve as a reference for researchers analyzing RDP1 and as a template for developing standardized QC workflows for other crop diversity panels.

## The Resource

### Dataset

The Rice Diversity Panel 1 consists of 413 *O. sativa* accessions genotyped with the 44K Affymetrix SNP array. The panel represents the five major varietal groups—indica (IND, n=87), temperate japonica (TEJ, n=96), tropical japonica (TRJ, n=97), aus (AUS, n=57), and aromatic (AROMATIC, n=14)—plus an admixed group (ADMIX, n=62). Accession passport data (NSFTV ID, GSOR ID, accession name, country of origin, subpopulation) are provided as part of the resource.

### Baseline QC Pipeline

The recommended QC pipeline, established through systematic evaluation (see Technical Validation), is:

```
Step 1: --mind 0.10          Remove samples with >10% missing genotypes
Step 2: --maf 0.05           Remove SNPs with minor allele frequency < 5%
Step 3: --geno 0.05          Remove SNPs missing in >5% of samples
Step 4: --indep-pairwise 50 5 0.2   LD pruning (50-SNP window, 5-SNP step, r² < 0.2)
```

This pipeline reduces the initial 413 accessions × 26,474 SNPs to **379 accessions × 1,187 SNPs**. The final dataset retains genome-wide coverage, represents all six subpopulations in proportions closely matching the original panel, and effectively eliminates background LD.

### Complete QC Scenario Dataset

To enable comparative analyses and allow researchers to evaluate the impact of alternative filtering decisions, the resource includes 16 systematically generated QC scenarios (Table 1). Each scenario varies one filtering parameter while holding all others at the baseline values, following a stepwise (not factorial) design:

| Dimension | Scenarios | Parameters varied |
|---|---|---|
| Sample filtering | 5 | mind = none, 0.10 (after QC), 0.10 (baseline), 0.05, 0.02 |
| MAF filtering | 5 | MAF = none, 0.01, 0.03, 0.05 (baseline), 0.10 |
| SNP missingness | 5 | geno = none, 0.20, 0.10, 0.05 (baseline), 0.02 |
| LD pruning | 5 | none, r² < 0.8, r² < 0.5, 50/5/0.2 (baseline), 100/10/0.2 |

For each scenario, the resource provides:
- PLINK binary files (.bed/.bim/.fam) for the filtered dataset
- ADMIXTURE output (K = 1–10, with Q matrices and CV error)
- PCA eigenvectors and eigenvalues (10 PCs)
- Allele frequency estimates (PLINK --freq)
- Hardy-Weinberg equilibrium statistics (PLINK --hardy)
- Pairwise LD estimates (PLINK --r2)
- FST values (PLINK --fst, Weir & Cockerham 1984)
- Individual-level genetic distance matrices (PLINK --distance)

### Population-Structure Inference Resource

For all 16 scenarios, ADMIXTURE was run at K = 1–10 with default convergence criteria. Cluster assignments at K = 5 (the biologically meaningful number of groups) are provided with label-aligned cross-references, enabling direct comparison of ancestry assignments across any pair of scenarios. PCA eigenvectors and DAPC results are also included.

### Excluded Accession Documentation

All 34 accessions removed during the baseline QC pipeline are documented with full passport metadata (Table S14). Four accessions were removed by sample-level missingness (mind > 0.10); thirty were removed during the stepwise SNP filtering pipeline as their missingness exceeded acceptable thresholds after SNP-level filters were applied. The distribution of excluded accessions across subpopulations is proportional to their representation in the original panel (Table 2).

## Technical Validation

### Dataset Composition After QC

The baseline pipeline removed 34 of 413 accessions (8.2%). Subpopulation representation before and after QC is summarized in Table 2:

| Subpopulation | 413 set | 379 set | Excluded | % change |
|---|---|---|---|---|
| ADMIX | 62 | 54 | 8 | −12.9% |
| AROMATIC | 14 | 13 | 1 | −7.1% |
| AUS | 57 | 49 | 8 | −14.0% |
| IND | 87 | 77 | 10 | −11.5% |
| TEJ | 96 | 93 | 3 | −3.1% |
| TRJ | 97 | 93 | 4 | −4.1% |

No subpopulation was disproportionately affected (χ² test, p > 0.05). Country representation was also preserved, with the same eight countries comprising the top sources of germplasm in both the 413- and 379-accession sets.

### Validation 1: Effects on SNP and Sample Retention

Sample filtering had the largest effect on dataset size (Figure 2). Mind = 0.02 retained only 98 accessions (24%) with 597 SNPs; mind = 0.05 retained 293 (71%) with 1,064 SNPs. Mind = 0.10 (recommended) retained 379 accessions (92%) with 1,187 SNPs.

LD pruning had the largest effect on SNP count. Without pruning, all 26,474 SNPs were retained. Pruning at r² < 0.8, r² < 0.5, and the recommended 50/5/0.2 window reduced SNP counts to 9,609, 4,278, and 1,187, respectively. MAF and GENO filters had moderate effects on SNP count (range: 844–1,638 SNPs) and negligible effects on sample retention.

**Key result**: The recommended thresholds retain 92% of accessions and 4.5% of SNPs, balancing data quality and genome coverage.

### Validation 2: Effects on Genetic Diversity

Mean PIC (polymorphism information content), expected heterozygosity (He), and the Shannon diversity index were nearly invariant across all 16 scenarios (PIC range: 0.27–0.31; He range: 0.27–0.31). The only notable deviation was MAF = 0.10, which reduced PIC to 0.27 by removing intermediate-frequency SNPs, and the absence of LD pruning, which slightly inflated PIC to 0.31.

Observed heterozygosity (Ho) was similarly stable (0.11–0.14) across most scenarios, with MAF = 0.10 producing an inflated value (0.19) due to the exclusion of low-MAF homozygous variants.

**Key result**: Diversity metrics per SNP are robust to QC variation. The primary effect of QC filtering is on the *number* of SNPs retained, not on their average diversity.

### Validation 3: Effects on Population Structure

**PCA stability.** PC1 and PC2 loadings from the baseline pipeline were highly correlated with those from all alternative scenarios (PC1 r > 0.99 for MAF, GENO, and LD scenarios). The only scenario with reduced PC1 correlation was the most stringent sample filter (mind = 0.02, 98 accessions; PC1 r = 0.82). PC2 correlations remained above 0.90 for all non-sample scenarios.

**ADMIXTURE cross-validation.** CV error curves exhibited the same qualitative shape across all 16 scenarios: rapid decrease from K = 1 to K = 2, plateau around K = 5–7, and a slight increase at K > 7. The CV minimum occurred at K = 5 for 10 of 16 scenarios, and at K = 4 for the remaining 6. The Evanno ΔK method peaked at K = 2 across all scenarios, reflecting the deepest hierarchical split in *O. sativa* (indica–japonica divergence). The K = 5 model is recommended as the biologically meaningful number of groups.

**Cluster stability.** After label alignment, mean cluster assignment agreement with the baseline was 93.1% across all 15 non-baseline scenarios (Figure 3). Agreement exceeded 99% for mind = none, mind = 0.10 after QC, and most MAF, GENO, and LD scenarios. The largest deviations occurred at extreme sample sizes: mind = 0.05 (293 accessions, 58.0% agreement) and mind = 0.02 (98 accessions, 56.1% agreement). Pairwise ARI and NMI values (Figures 4–5) confirmed this pattern, with MAF and GENO scenarios clustering together at ARI > 0.95 and sample filtering scenarios showing the greatest inter-scenario distance (ARI range: 0.10–0.58).

**Key result**: Population structure inference (PCA, ADMIXTURE) is robust across a wide range of QC thresholds when sample size exceeds ~350. Sample filtering below this threshold substantially alters ancestry assignments.

### Validation 4: Effects on FST and AMOVA

Mean FST across all 16 scenarios ranged from 0.29 to 0.61 (Figure 6). The most important validation result is the **stability of FST under recommended QC choices**: all scenarios with sample sizes >350 and standard MAF/GENO thresholds produced mean FST values between 0.39–0.45 (Table 3).

Notable deviations:
- **No LD pruning**: FST = 0.609 (+44% relative to baseline 0.423), demonstrating that LD pruning is essential for unbiased FST estimation
- **Mind = 0.02**: FST = 0.293 (−30.7%), driven by loss of rare or divergent accessions
- **MAF = none**: FST = 0.348 (−17.8%), as rare variants inflate within-population diversity

AMOVA results mirrored FST: among-population variance (ΦST) ranged from 16.1% (mind = 0.02) to 37.9% (no LD pruning), with the baseline at 30.2%.

**Key result**: FST from the recommended pipeline (0.423) is representative of estimates obtained across reasonable QC choices. LD pruning is the single most important factor for obtaining unbiased FST.

### Validation 5: Effects on LD Decay

LD decay was evaluated across five pruning strategies (Figure 7). Without pruning, the RDP1 dataset exhibited substantial LD (mean r² at 0–10 kb = 0.43; 48.1% of pairs with r² > 0.2), decaying to r² ≈ 0.1 at approximately 500–1,000 kb. This is consistent with the original RDP1 report (Zhao et al., 2011).

Progressive pruning reduced short-range LD:
| Pruning | SNPs | Mean r² (0–10 kb) | % pairs r² > 0.2 |
|---|---|---|---|
| None | 26,474 | 0.430 | 48.1% |
| r² < 0.8 | 9,609 | 0.207 | 23.0% |
| r² < 0.5 | 4,278 | 0.138 | 13.1% |
| **50/5/0.2 (recommended)** | **1,187** | **0.074** | **5.8%** |
| 100/10/0.2 | 822 | 0.070 | 1.4% |

The recommended 50/5/0.2 window reduces the SNP set by 95.5% while effectively eliminating background LD and retaining genome-wide coverage. At this threshold, r² ≈ 0.1 at approximately 388 kb, quantitatively consistent with the original RDP1 estimate of 500 kb–1 Mb (the difference attributable to SNP density, MAF threshold, and binning methodology).

**Key result**: The recommended pruning window effectively removes background LD. LD decay estimates from the curated dataset are reliable for GWAS marker spacing and population-genetic interpretation.

## Usage Notes

### How to Use This Resource

The curated dataset (379 accessions × 1,187 SNPs) and the complete 16-scenario dataset are available in PLINK binary format. Researchers can use the recommended baseline dataset directly or consult the sensitivity analysis to choose alternative thresholds appropriate for their specific research questions.

**For population-structure analysis**: Use the baseline dataset with ADMIXTURE K = 5. The Q matrices provided can be used directly for GWAS population-structure correction.

**For FST estimation**: LD pruning is essential. Use the baseline pruned dataset; unpruned data will substantially overestimate FST.

**For GWAS**: The baseline dataset is suitable for standard GWAS. Researchers requiring higher marker density can relax the LD pruning threshold to r² < 0.5 (4,278 SNPs) at the cost of slightly inflated FST (+2.6%).

**For diversity analysis**: The baseline dataset provides representative diversity estimates. MAF = 0.03 can be used to recover 17% more SNPs with negligible impact on diversity or FST (−6.9%).

### Applying the Pipeline to New Data

The stepwise QC pipeline and comparison framework can be applied to any SNP dataset. The R scripts provided in the repository include all steps: PLINK-based filtering, ADMIXTURE batch execution, PCA, FST, AMOVA, genetic diversity calculation, LD decay analysis, and sensitivity metrics (confusion matrices, ARI/NMI, PC correlation). The workflow diagram (Figure 1) provides a visual guide to the pipeline structure.

### Comparison with the Original RDP1

Users comparing results with the original RDP1 publication should note that:
1. The baseline dataset (379 accessions) is 8% smaller than the original 413-accession panel
2. The excluded accessions are documented in Table S14
3. LD decay estimates are quantitatively similar (r² ≈ 0.1 at ~388 kb vs. 500 kb–1 Mb in the original)
4. Population structure at K = 5 is consistent with the original five-group classification

## Data Availability

All data and code are available at:
- **GitHub repository**: https://github.com/yehiakhidr-byte/rice
- **Zenodo archive**: [DOI to be added upon publication]

The repository contains:
- 40+ R and shell scripts for the complete analysis pipeline
- PLINK binary files for all 16 QC scenarios
- ADMIXTURE output (K = 1–10) for all scenarios
- PCA eigenvectors and eigenvalues
- Allele frequency and HWE statistics
- FST and AMOVA results
- LD decay data
- All supplementary figures and tables
- Full RDP1 passport metadata and cross-reference tables
- Complete parameter documentation
- Workflow diagram

All analyses were performed with publicly available software: PLINK 1.90b7.2, ADMIXTURE v1.3.0, and R 4.x (packages: adegenet, pegas, poppr, vegan, ggplot2).

## Conclusions

RiceDP1-QC provides the first systematically validated, fully reproducible QC pipeline for the Rice Diversity Panel 1. The resource addresses a long-standing gap in RDP1 research by establishing empirically grounded threshold recommendations, quantifying the sensitivity of key population-genetic metrics to QC decisions, and providing a complete computational workflow that can be directly applied or adapted by other researchers.

The principal findings of the technical validation are:
1. **The recommended pipeline (mind 0.10 → MAF 0.05 → geno 0.05 → LD 50/5/0.2) produces robust, reproducible population-genetic inferences.** FST, diversity, and cluster assignments are stable across a wide range of alternative thresholds.
2. **Sample size is the most critical variable.** Datasets with <350 accessions produce substantially different FST and ancestry assignments.
3. **LD pruning is essential for unbiased FST estimation.** Unpruned data inflates FST by 44%.
4. **All methods (PCA, ADMIXTURE, DAPC) converge on the same five-group structure** when sample size is adequate.

By providing a fully documented, version-controlled, and citable resource, RiceDP1-QC aims to improve the reproducibility and comparability of RDP1-based analyses and to serve as a template for similar efforts in other crop diversity panels.