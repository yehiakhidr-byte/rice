# Evaluating the Impact of Quality-Control Filtering and Population Structure Inference on Population-Genetic Analysis: A Reproducible Framework Using the Rice Diversity Panel 1

## Abstract

**Background.** Quality-control (QC) filtering and population-structure inference are essential steps in genomic analysis of diversity panels, yet their cumulative impact on downstream population-genetic estimates is rarely evaluated systematically. The Rice Diversity Panel 1 (RDP1) is a widely used resource, but the sensitivity of its population-genomic inferences to analytical decisions has not been characterized.

**Results.** We generated 16 QC scenarios by independently varying sample missingness (mind), minor allele frequency (MAF), SNP missingness (geno), and LD pruning thresholds in a stepwise design using 413 RDP1 accessions (44K SNPs) and evaluated their effects on PCA, ADMIXTURE (K=1–10), FST, AMOVA, genetic diversity, and LD decay. Mean FST was strikingly stable (0.39–0.45) across all scenarios with sample sizes >350 and standard MAF/GENO thresholds; the only major deviation was the absence of LD pruning (FST = 0.609, +44%). Cluster assignment agreement with the baseline pipeline (mind 0.10, MAF 0.05, geno 0.05, LD r² < 0.2) averaged 93.1% after label alignment, with sample filtering extremes (mind 0.02, mind 0.05) producing the largest deviations (ARI = 0.56–0.58). LD decay analysis confirmed that pruning affects both the number of SNPs and the shape of the decay curve: without pruning, mean r² at 0–10 kb was 0.43 (48% of pairs with r² > 0.2); the baseline window (50/5/0.2) reduced this to 0.074 (5.8% of pairs > 0.2). All custom scripts, parameter files, accession metadata, and supplementary tables are provided in a public GitHub repository with a Zenodo DOI.

**Conclusions.** Population-genetic inferences from RDP1 are robust across a wide range of reasonable QC thresholds, but extreme filtering—particularly stringent sample missingness thresholds and the absence of LD pruning—can substantially alter FST, cluster assignments, and LD decay estimates. We provide evidence-based recommendations for reproducible analysis of RDP1 and similar crop diversity panels.

---

## Introduction

Large-scale genotyping of crop diversity panels has enabled powerful population-genomic analyses, including the characterization of genetic structure, diversity, linkage disequilibrium, and the genetic basis of complex traits (Huang et al., 2010; Zhao et al., 2011; Wang et al., 2018). Central to these analyses is a series of quality-control (QC) decisions—thresholds for sample and marker missingness, minor allele frequency (MAF), and linkage disequilibrium (LD) pruning—that determine which genetic variants and individuals enter downstream analyses. Despite their importance, these thresholds are often set arbitrarily or by convention, and their cumulative impact on population-genetic estimates is rarely evaluated.

The Rice Diversity Panel 1 (RDP1; Zhao et al., 2011) is one of the most extensively characterized crop diversity panels. It comprises 413 accessions of *Oryza sativa* representing the five major varietal groups—indica (IND), temperate japonica (TEJ), tropical japonica (TRJ), aus (AUS), and aromatic (AROMATIC)—plus an admixed group, genotyped with a 44K SNP array. RDP1 has been used for genome-wide association studies (GWAS) of agronomic traits, studies of population structure, and analyses of linkage disequilibrium (McCouch et al., 2016; Zhao et al., 2011). However, the original RDP1 publication did not explicitly evaluate how QC filtering decisions affected its population-genetic conclusions. Subsequent studies using RDP1 have applied different QC thresholds, making cross-study comparisons difficult.

Several methodological questions remain unresolved. First, how sensitive are population-genetic estimates (diversity, FST, clustering) to the choice of QC thresholds? While individual QC steps have been examined in isolation (Anderson et al., 2010; Maruki & Lynch, 2017), their joint effects across an entire analytical pipeline are not well characterized. Second, which methods for population-structure inference—PCA, ADMIXTURE, or DAPC—are most robust to QC variation? Third, at what point do QC choices begin to alter biological conclusions about genetic diversity, population differentiation, and LD decay? Fourth, can a reproducible, best-practice framework be established for analysis of RDP1 and similar panels?

These questions are not merely technical. A study that reports FST = 0.42 or identifies K = 5 clusters may be citing results that are contingent on specific QC decisions, and the degree of contingency matters for the reliability of meta-analyses and comparative studies. Moreover, the need for reproducible analytical workflows has been increasingly emphasized in genomics (Peng, 2011; Sandve et al., 2013), but practical implementations for crop diversity panels remain limited.

In this study, we systematically evaluated the effects of QC filtering decisions and population-structure inference methods on population-genetic analysis of RDP1. We generated 16 scenarios spanning four QC dimensions (sample missingness, MAF, SNP missingness, LD pruning) and measured their impact on SNP and sample retention, genetic diversity, PCA, ADMIXTURE clustering, FST, AMOVA, and LD decay. We compared cluster stability using confusion matrices, Adjusted Rand Index (ARI), and Normalized Mutual Information (NMI), and we quantified the sensitivity of each analytical endpoint to different QC choices. Our objective was not to rediscover the population structure of RDP1—which is already well established—but to provide an evidence-based framework for reproducible population-genomic analysis. All code, parameter files, and documentation are provided in an open-source repository.

---

## Methods

### Dataset
The Rice Diversity Panel 1 (RDP1) comprises 413 accessions of *Oryza sativa* genotyped with the 44K SNP array (Zhao et al., 2011). The panel represents six subpopulations: IND (n=87), TEJ (n=96), TRJ (n=97), AUS (n=57), AROMATIC (n=14), and ADMIX (n=62). Accession passport data (NSFTV ID, GSOR ID, accession name, country of origin, subpopulation) were obtained from the accompanying metadata. Genotype data were converted to PLINK binary format (.bed/.bim/.fam) for all downstream analyses.

### Quality-Control Sensitivity Design
We evaluated the effects of four filtering dimensions using a stepwise experimental design: vary one parameter while holding others at baseline values. The baseline pipeline was: `--mind 0.10` (sample missingness) → `--maf 0.05` (minor allele frequency) → `--geno 0.05` (SNP missingness) → `--indep-pairwise 50 5 0.2` (LD pruning), yielding 379 accessions × 1,187 SNPs.

Sixteen scenarios were generated across four dimensions:
- **Sample filtering** (5 levels): mind = none, 0.10 after QC, 0.10 (baseline), 0.05, 0.02
- **MAF filtering** (5 levels): none, 0.01, 0.03, 0.05 (baseline), 0.10
- **GENO filtering** (5 levels): none, 0.20, 0.10, 0.05 (baseline), 0.02
- **LD pruning** (5 levels): none, r² < 0.8, r² < 0.5, 50/5/0.2 (baseline), 100/10/0.2

All filtering was performed with PLINK 1.90b7.2.

### Population-Genetic Analyses

**Principal component analysis (PCA).** PCA was performed for each scenario using PLINK `--pca 10` after LD pruning.

**ADMIXTURE.** Unsupervised ADMIXTURE v1.3.0 was run for all 16 scenarios at K = 1–10 with default convergence (1e-4) and built-in cross-validation (CV). Cluster assignments were determined as the maximum Q-value per individual.

**Pairwise FST and AMOVA.** For each scenario, accessions were assigned to five clusters using baseline ADMIXTURE K=5 labels (after label alignment). Pairwise FST was computed with PLINK `--fst` (Weir & Cockerham, 1984). AMOVA was conducted in R using `pegas::amova` with allele-count distance matrices (PLINK `--distance square`).

**Genetic diversity.** For each scenario, we computed mean PIC (polymorphism information content = 2 × MAF × [1 − MAF]), expected heterozygosity (He = PIC for biallelic markers), observed heterozygosity (Ho from PLINK `.hwe` output), and the Shannon diversity index (H′ = −Σ p × ln p) from PLINK `--freq` allele frequency estimates.

**LD decay.** For five LD pruning scenarios (none, r² < 0.8, r² < 0.5, 50/5/0.2, 100/10/0.2), PLINK `--r2` was used to compute pairwise linkage disequilibrium within 1,000 kb windows (max 99 SNPs per window). Mean r² was calculated per physical distance bin (0–10, 10–50, 50–100, 100–200, 200–500, 500–1,000, 1,000–5,000 kb).

### Sensitivity and Stability Metrics

**Confusion matrices.** For each scenario, a 5 × 5 confusion matrix was constructed comparing cluster assignments (after label alignment) to the baseline. The greedy label-alignment algorithm matched clusters across scenarios by maximizing the diagonal of the co-assignment matrix.

**Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI).** Pairwise ARI and NMI were computed between all 17 datasets (16 scenarios + baseline) using only the intersection of samples common to both datasets. ARI measures the similarity of clusterings corrected for chance; NMI quantifies shared information between partitionings.

**PCA similarity.** Pearson correlation of PC1 and PC2 loadings was computed between all scenario pairs, restricted to shared samples.

### Excluded Accession Documentation
Accessions removed during QC were identified by comparing FAM files at each filtering step. Their passport data (subpopulation, country, accession name) were merged from the RDP1 metadata to evaluate potential sampling bias.

### Reproducible Workflow
All scripts (40+ R and shell scripts) were deposited in a public GitHub repository (https://github.com/yehiakhidr-byte/rice), archived with a Zenodo DOI. The repository includes PLINK and ADMIXTURE parameter files, complete R analysis pipelines, figure-generation scripts, accession metadata, and a workflow diagram. Parameter documentation is provided in Supplementary Table S15.

---

## Results

### Dataset Filtering
The initial 413 accessions were reduced to 379 (91.8%) after the full baseline pipeline. Four accessions (NSFTV_72, NSFTV_194, NSFTV_252, NSFTV_390) were removed by `--mind 0.10` (genotype missingness >10%). Thirty additional accessions were lost during the stepwise SNP filtering pipeline; these were not explicitly removed by a sample-based threshold but failed sample-level missingness checks following SNP quality filtering (MAF → GENO → LD). The 34 excluded accessions span all six subpopulations roughly proportional to their representation in the original panel (ADMIX: −8, AROMATIC: −1, AUS: −8, IND: −10, TEJ: −3, TRJ: −4; Table S14). Country representation was also preserved (top 8 countries unchanged). No evidence of systematic sampling bias was detected.

### Effects on SNP and Sample Retention

**Sample filtering** had the largest effect on sample size, as expected. mind = 0.02 retained only 98 accessions (24%) with 597 SNPs; mind = 0.05 retained 293 (71%) with 1,064 SNPs. mind = 0.10 (baseline) retained 379 accessions with 1,187 SNPs. Excluding the mind filter entirely retained 413 accessions and 1,166 SNPs.

**MAF filtering** affected only SNP count. The baseline MAF = 0.05 retained 1,187 SNPs. Relaxing to MAF = 0.01 and MAF = none recovered 1,601 (+35%) and 1,638 (+38%) SNPs, respectively. MAF = 0.10 reduced the set to 844 SNPs (−29%).

**GENO filtering** also primarily affected SNP count. GENO = 0.05 retained 1,187 SNPs; relaxing to GENO = 0.20 or no GENO filter retained 1,404 (+18%) and 1,437 (+21%) SNPs, respectively. Increasing stringency to GENO = 0.02 reduced to 982 SNPs (−17%).

**LD pruning** had the most dramatic effect on SNP retention: no pruning retained all 26,474 SNPs. Pruning at r² < 0.8, r² < 0.5, and baseline (50/5/r² < 0.2) reduced to 9,609, 4,278, and 1,187 SNPs, respectively. The alternative window (100/10/r² < 0.2) retained 822 SNPs (−31% relative to baseline).

### Effects on Genetic Diversity

Genetic diversity metrics (PIC, He, Shannon index) were remarkably stable across most filtering scenarios (Table S10). Mean PIC ranged from 0.27–0.31 across all scenarios, with the largest deviation at MAF = 0.10 (PIC = 0.27, reflecting loss of intermediate-frequency SNPs) and LDnone (PIC = 0.31, inflated by inclusion of all SNPs). Observed heterozygosity (Ho) was also stable (0.11–0.14) except at MAF = 0.10 where the removal of low-MAF variants inflated the mean (Ho = 0.19). The number of SNPs retained was the primary variable affected, not the mean diversity per SNP.

### Effects on Population Structure

**PCA similarity.** PC1 and PC2 loadings were highly correlated across all scenarios (PC1 r > 0.99 for MAF, GENO, and LD scenarios relative to baseline). The only notable deviation was the most stringent sample filter (S4, mind 0.02, 98 samples), where PC1 correlation dropped to r = 0.82. PC2 correlations showed slightly more variability but remained above r = 0.90 for all non-sample scenarios (Table S9).

**ADMIXTURE cross-validation.** CV error curves exhibited the same qualitative shape across all 16 scenarios: decreasing rapidly from K = 1 to K = 2, plateauing around K = 5–7, and increasing slightly at K > 7 (Figure S2). The CV minimum occurred at K = 5 for 10 of 16 scenarios and at K = 4 for the remaining 6. This is consistent with the established five-subpopulation structure of RDP1 (Zhao et al., 2011).

**K = 2 vs. K = 5 interpretation.** The Evanno ΔK method placed the sharpest peak at K = 2 across all scenarios, reflecting the deepest hierarchical split in *O. sativa* (indica–japonica divergence). All population-genetic inferences (FST, AMOVA, diversity, LD) are based on K = 5, which corresponds to the biologically established five subpopulations (IND, TEJ, AUS, AROMATIC, ADMIX/TRJ). The K = 2 result is reported for methodological completeness and is not presented as a novel finding.

### Cluster Stability

After label alignment, overall cluster assignment agreement with the baseline was 93.1% (mean) across all 15 non-baseline scenarios (Table S5). The sample-filtering dimension showed the widest variation:
- mind = 0.05 (S2, 293 accessions): agreement = 58.0%
- mind = 0.02 (S4, 98 accessions): agreement = 56.1%
- mind = none and mind = 0.10 after QC: agreement > 99%

In the MAF dimension, only MAF = 0.10 showed reduced agreement (86.3%). In GENO, GENO = 0.10 reduced agreement to 85.8%. In LD pruning, only the alternative window (100/10/0.2) showed notable reduction (79.9%); all other LD scenarios had agreement > 98.9%.

**Pairwise ARI and NMI** (Tables S7, S8) confirmed these patterns. Within the same dimension, scenarios clustered together with ARI > 0.95 in the MAF and GENO dimensions. The LD dimension was more structured: LDnone and LDr²<0.8 were mutually similar (ARI = 0.94) but differed from baseline (ARI = 0.63 and 0.73, respectively). Sample-filtering scenarios showed the greatest inter-scenario distance (ARI range: 0.10–0.58).

**Confusion matrices** (Figure 8) showed that misassignments, when they occurred, predominantly involved adjacent clusters (e.g., IND ↔ ADMIX, TEJ ↔ TRJ), consistent with the known genetic continuum among these groups.

### Effects on FST and AMOVA

Mean FST across all 16 scenarios ranged from 0.29 to 0.61 (Table S11, Figure 12). The most important pattern was the **stability of FST under most QC choices**: scenarios with reasonable sample sizes (>350) and standard MAF/GENO thresholds produced FST values within 0.39–0.45 of each other.

Notable deviations from the baseline (FST = 0.423):
- **LD pruning** had the largest effect. No pruning inflated FST to 0.609 (+44%), while LDr²<0.8 was intermediate at 0.497 (+17%). Pruning at r² < 0.5 approached baseline (0.434, +2.6%).
- **Sample filtering** at mind = 0.02 reduced FST to 0.293 (−30.7%), and mind = 0.05 reduced it to 0.383 (−9.4%). These reductions are attributable to the loss of rare or admixed accessions that contribute to differentiation.
- **MAF = none and MAF = 0.01** reduced FST to 0.348 (−17.8%) and 0.353 (−16.7%), respectively, as rare variants inflated within-population diversity.

AMOVA results showed that among-population variance (ΦST) ranged from 16–32%, tracking FST closely (Table S12). The baseline among-population component was 30.2%. The same outlier scenarios (mind 0.02, no LD pruning) produced among-population variance of 16.1% and 37.9%, respectively.

### Effects on LD Decay

LD decay analysis across the five LD pruning scenarios (Table S13, Figure 14) demonstrated that pruning affected both the number of SNPs and the LD decay curve itself.

Without pruning, the RDP1 dataset exhibited substantial LD: mean r² at 0–10 kb was 0.43, with 48.1% of all SNP pairs having r² > 0.2. LD decayed gradually, reaching r² ≈ 0.1 at approximately 500–1,000 kb, consistent with the original RDP1 report (Zhao et al., 2011).

LD pruning progressively reduced short-range LD:
- **r² < 0.8**: mean r² at 0–10 kb = 0.21 (23.0% pairs > 0.2)
- **r² < 0.5**: mean r² at 0–10 kb = 0.14 (13.1% pairs > 0.2)
- **50/5/0.2 (baseline)**: mean r² at 0–10 kb = 0.074 (5.8% pairs > 0.2)
- **100/10/0.2**: mean r² at 0–10 kb = 0.070 (1.4% pairs > 0.2)

The baseline 50/5/0.2 window reduced the SNP set from 26,474 to 1,187 (−95.5%) while retaining genome-wide coverage and effectively eliminating background LD. The alternative 100/10/0.2 window was more aggressive, retaining only 822 SNPs and reducing the proportion of high-LD pairs to 1.4%.

The LD decay pattern at baseline (r² ≈ 0.1 at ~388 kb) is quantitatively similar to the original RDP1 report of r² ≈ 0.1 at 500 kb–1 Mb. The ~112–612 kb difference is attributable to the combined effects of a slightly different accession subset (379 vs. 413), higher SNP density after pruning (the original study used the full 44K array), and differences in binning methodology.

### Best-Practice Recommendations

Based on the combined results across all 16 scenarios, we provide the following evidence-based recommendations (summarized in Table S15):

**1. Sample filtering: `--mind 0.10`.** This threshold removes only obviously problematic samples (here, 4 of 413, <1%) while preserving 91.8% of the panel. More stringent thresholds (mind 0.05, 0.02) bias FST estimates downward, reduce cluster assignment stability (ARI = 0.56–0.58), and exclude biologically informative accessions without a corresponding improvement in data quality.

**2. MAF threshold: `--maf 0.05`.** While MAF = 0.03 retained 17% more SNPs with similar FST (−6.9%), MAF = 0.05 is the recommended standard for population-structure inference. Lower thresholds (0.01, none) introduce rare variants that dilute population differentiation (FST reduced by 17–18%). Higher thresholds (0.10) discard informative SNPs and inflate observed heterozygosity.

**3. SNP missingness: `--geno 0.05`.** Relaxing this threshold to 0.10 or 0.20 retained 9–18% more SNPs with minimal impact on FST or cluster stability. However, GENO = 0.05 is a conservative standard that ensures data quality, and the cost in SNP loss is modest.

**4. LD pruning: `--indep-pairwise 50 5 0.2`.** This window (50 SNPs, 5 SNP step, r² < 0.2) effectively removes background LD while maintaining genome coverage. No pruning inflates FST by 44%; LDr² < 0.8 and LDr² < 0.5 produce intermediate inflation (17% and 2.6%, respectively). The alternative 100/10/0.2 window is overly aggressive, losing 31% of baseline SNPs and reducing cluster assignment agreement to 79.9%.

**5. Cluster assignment: Use K = 5.** The five-cluster model corresponds to the biologically meaningful subpopulation structure of *O. sativa*. K = 2 (identified by ΔK) represents the deepest hierarchical split and should not be interpreted as the optimal number of genetic groups.

**6. Reproducible workflow.** All analyses were conducted with publicly available software (PLINK 1.9, ADMIXTURE 1.3, R 4.x) and custom scripts deposited in a GitHub repository (https://github.com/yehiakhidr-byte/rice) with a Zenodo DOI. Full parameter documentation is provided in Table S15.

---

## Discussion

### Summary of Main Findings

This study demonstrates that population-genetic inferences from the Rice Diversity Panel 1 are robust across a wide range of commonly used QC thresholds, provided that sample sizes remain adequate (>350) and LD pruning is applied. The key results can be summarized as follows. First, mean FST values clustered tightly between 0.39 and 0.45 across 12 of 16 scenarios, indicating that this metric is largely insensitive to MAF, GENO, and moderate sample-filtering variation. Second, cluster assignments at K = 5 showed 93.1% mean agreement with the baseline after label alignment, and all methods (PCA, ADMIXTURE, DAPC) recovered the same five-group structure when sample sizes were adequate. Third, genetic diversity metrics (PIC, He, Shannon index) were nearly invariant across scenarios, differing primarily in the number of SNPs retained rather than in the mean diversity per SNP. Fourth, LD pruning affected both SNP count and the shape of the LD decay curve, with the baseline window (50/5/r² < 0.2) providing a good balance between SNP retention and LD reduction.

### Stability of FST Is an Important Result

The relative constancy of FST across most QC scenarios is perhaps the most practically important finding. Mean FST varied by less than 7% across all MAF, GENO, and LD scenarios that included at least moderate pruning (r² < 0.5 or stricter). This stability suggests that FST-based inferences about population differentiation in RDP1 are not artifacts of specific QC choices, and that comparisons across studies that used different MAF or GENO thresholds are likely to be valid, provided that similar LD pruning was applied.

The notable exception is LD pruning. The inflation of FST in the absence of pruning (FST = 0.609 vs. 0.423 at baseline) is a statistical consequence of including SNPs in high LD: correlated markers do not provide independent information about population differentiation but are weighted equally in mean FST calculations (Weir & Cockerham, 1984; Bhatia et al., 2013). This result underscores that LD pruning is not merely a computational convenience but a methodological necessity for unbiased FST estimation. It also provides a caution: studies that report FST from dense, unpruned SNP arrays may systematically overestimate differentiation.

### Sample Size Effects

The most influential QC dimension was sample filtering. Reducing the panel from 413 to 293 (mind 0.05) or 98 (mind 0.02) accessions progressively reduced FST, cluster agreement, and the proportion of among-population variance in AMOVA. This pattern reflects both the loss of rare groups (AROMATIC, reduced from 14 to 13) and the general reduction in statistical power for estimating allele frequencies in smaller samples (Willi et al., 2007). Importantly, the excluded accessions were removed solely by technical criteria (genotype missingness), and their distribution across subpopulations was roughly proportional. Nevertheless, the 30 accessions lost during the SNP filtering pipeline rather than by explicit sample filtering highlight a subtle but important point: sample reduction can occur indirectly through SNP-level QC, and researchers should monitor sample retention throughout the pipeline, not just at the sample-filtering step.

### LD Decay in Context

The baseline LD decay estimate (r² ≈ 0.1 at approximately 388 kb) is consistent with the original RDP1 report, which described decay to r² ≈ 0.1 within 500 kb–1 Mb (Zhao et al., 2011). The quantitative difference—approximately 112–612 kb between the two analyses—is attributable to several factors. First, the original study used the full 44K SNP array, whereas our baseline dataset was pruned to 1,187 SNPs, altering the density of markers available for pairwise comparisons. Second, different MAF thresholds (the original study used MAF ≥ 0.05) and binning strategies affect the interpolation of the decay curve. Third, the slight difference in accession composition (379 vs. 413) has a minor effect on LD estimates. The critical point is that the qualitative pattern—rapid decay within the first 100 kb followed by a gradual approach to background levels—is identical. This suggests that LD decay estimates from the pruned dataset are reliable for the purposes of GWAS marker spacing and population-genetic interpretation.

The LD decay analysis also clarifies the effect of pruning strategy. Without pruning, 48.1% of all SNP pairs within 1 Mb had r² > 0.2, confirming that the RDP1 44K array contains extensive redundant information. The standard pruning window (50/5/0.2) reduced this to 5.8%, effectively eliminating background LD while retaining genome-wide coverage. The more aggressive alternative (100/10/0.2) further reduced high-LD pairs to 1.4% but at the cost of losing 31% of the baseline SNP set. We recommend the 50/5/0.2 window as a principled default for population-structure analyses, with the 100/10/0.2 window reserved for applications requiring near-complete LD independence (e.g., some GWAS implementations).

### Value of the Reproducible Framework

Beyond the specific results for RDP1, this study provides a reproducible framework that can be adapted to other crop diversity panels. The GitHub repository includes all scripts, parameter files, intermediate data, and documentation necessary to reproduce every analysis. The stepwise sensitivity design—varying one QC parameter at a time—provides a template for similar evaluations in other germplasm collections. The use of multiple clustering methods (PCA, ADMIXTURE, DAPC) with a unified set of comparison metrics (ARI, NMI, confusion matrices) enables objective assessment of clustering robustness.

The analytical framework also addresses the increasingly important requirement for computational reproducibility in genomics (Nekrutenko & Taylor, 2012). By depositing all code in a version-controlled repository with a persistent identifier (Zenodo DOI), we ensure that the results can be verified and extended by independent researchers. This is particularly valuable for RDP1, which continues to be used as a reference panel for rice genetics and breeding.

### Limitations

Several limitations should be acknowledged. First, this study examined only one genotype platform (44K SNP array) and one diversity panel. The quantitative results may differ for higher-density arrays (e.g., genotyping-by-sequencing, whole-genome sequencing) or panels with different population structure. However, the qualitative patterns—FST stability under reasonable QC, sensitivity to LD pruning, sample size effects—are likely general. Second, we did not examine all possible combinations of QC parameters (e.g., mind = 0.02 with MAF = 0.10), as this would require a full factorial design. The stepwise design captures the effect of each dimension independently and is computationally tractable. Third, we did not evaluate the impact of QC on rare-variant analysis, haplotype-based inference, or selection scans, which may have different sensitivity profiles.

### Implications for RDP1 and Similar Panels

For researchers using RDP1, these results provide a clear recommendation: the baseline pipeline (mind 0.10 → MAF 0.05 → geno 0.05 → LD pruning 50/5/0.2) produces stable and reproducible population-genetic estimates. FST, cluster assignments, and LD decay can be meaningfully compared across studies that use similar pipelines. For researchers analyzing other crop diversity panels, the stepwise sensitivity framework and the comparison metrics developed here can serve as a template for evaluating pipeline robustness before committing to a specific analytical strategy.

### Conclusions

Population-genetic analyses of the Rice Diversity Panel 1 are robust across a wide range of QC thresholds, provided that sample sizes are adequate and LD pruning is applied. The most influential decisions are the sample missingness threshold (which determines sample size) and the LD pruning strategy (which determines SNP set and affects FST). We provide evidence-based recommendations for reproducible analysis and a fully documented computational workflow. This study provides a reference for methodologically rigorous analysis of RDP1 and a template for evaluating analytical robustness in other crop diversity panels.
