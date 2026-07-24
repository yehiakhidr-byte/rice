## Recommended Figures and Tables for the Paper

### Main Figures (6)
| # | File | Title | Why in Main |
|---|------|-------|-------------|
| 1 | `Figure_1_workflow.pdf` | Experimental workflow and reproducible framework | Essential for the reviewer's reproducibility concern. Shows the stepwise QC design (16 scenarios × 4 dimensions) → population analyses → sensitivity comparison → recommendations. Without this, the reader cannot follow the paper's structure. |
| 2 | `Figure_6_QC_impact.pdf` | Impact of QC filtering on dataset composition | First-order effect: how many SNPs and samples survive each threshold. The reader needs to see this before any downstream results. |
| 3 | `Figure_8_confusion_matrix.pdf` | Cluster assignment consistency across QC scenarios | Directly addresses "did accessions remain indica/japonica/aus/aromatic/admixed or change?" Shows diagonal = stable, off-diagonal = switches. This is the strongest visual of cluster stability. |
| 4 | `Figure_12_FST_comparison.pdf` | Effects of QC filtering on mean FST | Most important quantitative result. The stability of FST (0.39–0.45) across reasonable QC is a publishable finding. The LDnone inflation (+44%) is a clear caution. This figure alone justifies the paper. |
| 5 | `Figure_14_LD_decay.pdf` | LD decay under different pruning strategies | Directly answers "does pruning change LD decay or only SNP count?" Shows both effects. Needed to contextualize against RDP1 (500 kb–1 Mb). |
| 6 | `Figure_16_recommendations.pdf` | Best-practice recommendations dashboard | The actionable output of the paper. Four panels show SNP retention, FST, ARI, and LD decay across all thresholds. Transforms the paper from "we compared methods" into "here is what to use." |

### Supplementary Figures (5)
| # | File | Title | Why Supplementary |
|---|------|-------|------------------|
| S1 | `Figure_PCA_structure_validation.pdf` | PCA of RDP1 colored by ADMIXTURE K=5 | Background — shows the well-known RDP1 structure. Readers familiar with RDP1 do not need this, but new readers do. |
| S2 | `Figure_S2_sensitivity_CV_curves.pdf` | ADMIXTURE CV error curves (all 16 scenarios) | Supporting evidence that K=5 is robust across QC choices. Too many panels for main text. |
| S3 | `Figure_9_ARI_NMI.pdf` | Pairwise ARI and NMI heatmaps | Comprehensive pairwise comparison. The key information (range of agreement) can be summarized in one sentence in main text. The full matrix belongs in supplementary. |
| S4 | `Figure_10_PCA_correlations.pdf` | PC1/PC2 correlation across scenarios | Supporting — shows PCA is robust. A single sentence ("PC1 r > 0.99 for all scenarios except mind 0.02") suffices in main text. |
| S5 | `Figure_7_cluster_stability.pdf` | Distribution of per-sample agreement scores | Supports the cluster stability narrative. The mean (93.1%) is reported in main text; the full distribution is supplementary. |

### Excluded Figures (not recommended)
| File | Why Exclude |
|------|-------------|
| `Figure_5_workflow.pdf` | Superseded by `Figure_1_workflow.pdf` (higher quality version) |
| `Figure_10_diversity.pdf` | Older version, superseded by `Figure_11_genetic_diversity.pdf` |
| `Figure_11_FST_comparison.pdf` | Older version, superseded by the newer `Figure_12_FST_comparison.pdf` |
| `Figure_13_AMOVA.pdf` | AMOVA tracks FST closely (r² > 0.99). A single sentence in main text ("AMOVA among-population variance ranged from 16–32%, closely tracking FST") suffices. |
| `Figure_S4_sensitivity_summary.pdf` | Redundant — the recommendation dashboard (Fig 16) conveys the summary more effectively |
| `Figure_S5_cluster_flow.pdf` | Redundant with confusion matrices (Fig 3) — shows the same data in a different format |

### Main Tables (2)
| Table | File | Why in Main |
|-------|------|-------------|
| 1 | `Table_S14_excluded_accessions.csv` | Directly addresses the reviewer's criticism about undocumented sample reduction. Must be in main text to demonstrate no sampling bias. |
| 2 | `Table_S15_recommendations.csv` | The actionable output. Summarizes the recommended thresholds with empirical impact on all metrics. |

### Supplementary Tables (13)
All others (S1–S13) belong in supplementary as they contain the full data underlying the figures. They are essential for reproducibility but would overwhelm the main text.

---

## Checks Against Paper Structure

| Paper section | Main Figures | Main Tables | Supplementary |
|---------------|-------------|-------------|---------------|
| Introduction — RDP1 background | — | — | Fig. S1 (PCA) |
| Methods — QC design | Fig. 1 (workflow) | — | — |
| Methods — population analyses | Fig. 1 (workflow) | — | — |
| Results — dataset filtering | Fig. 2 (QC impact) | Table 1 (excluded accessions) | Tables S1–S4 |
| Results — cluster stability | Fig. 3 (confusion matrices) | — | Figs. S3, S5; Tables S5–S8 |
| Results — FST and AMOVA | Fig. 4 (FST) | — | Fig. S3; Tables S11–S12 |
| Results — genetic diversity | (1–2 sentences in text) | — | Fig. S3; Table S10 |
| Results — PCA similarity | (1 sentence in text) | — | Fig. S4; Table S9 |
| Results — LD decay | Fig. 5 (LD decay) | — | Table S13 |
| Discussion — recommendations | Fig. 6 (dashboard) | Table 2 (recommendations) | — |
| Methods — reproducibility | Fig. 1 (GitHub/Zenodo) | — | All tables |

---

## Additional Checks

### What is missing?
1. **No phenotypic data for GWAS** — the paper outline mentions GWAS as a downstream component, but no phenotype data is available in the current dataset. Without it, the "consequences for GWAS" section cannot be written. Consider either: (a) obtaining phenotype data (e.g., flowering time, grain size from McCouch et al. 2016), or (b) dropping the GWAS framing and keeping the paper focused on population-genetic inference.

2. **No runtime benchmark table** — reviewers may ask "how long did ADMIXTURE take for each scenario?" This is easy to generate from the log files.

3. **No comparison with STRUCTURE** — the paper outline mentions STRUCTURE but we only ran ADMIXTURE. ADMIXTURE is computationally equivalent to STRUCTURE with the no-linkage model and is the standard for SNP-scale data. This should be explicitly stated in Methods to avoid reviewer questions.

### What is redundant?
- `Figure_13_AMOVA.pdf` adds no new information beyond FST
- `Figure_S4_sensitivity_summary.pdf` duplicates the recommendation dashboard
- `Table_S3_PCA_procrustes.csv` is superseded by the PC correlation table
- `Table_S6_ARI_matrix.csv` is identical information to `Table_S7_pairwise_ARI.csv` (different sorting)
- `Table_S11b_LD_snp_counts.csv` is subsumed by Table S13

### File cleanup recommendation
Delete these legacy/duplicate files from the repository to avoid confusion:
- `Figure_5_workflow.pdf`, `Figure_10_diversity.pdf`, `Figure_11_FST_comparison.pdf`, `Figure_13_AMOVA.pdf`
- `Figure_8_ARI_heatmap.pdf`, `Figure_9_NMI_heatmap.pdf` (superseded by combined `Figure_9_ARI_NMI.pdf`)
- `Table_S6_ARI_matrix.csv`, `Table_S7_NMI_matrix.csv`, `Table_S8_PCA_correlation.csv`, `Table_S10_AMOVA.csv`, `Table_S11_LD_decay.csv`, `Table_S11b_LD_snp_counts.csv`