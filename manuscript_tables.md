## Tables

### Table 1. Dataset Characteristics After Baseline QC
| Metric | 44K | HDRA |
|--------|-----|------|
| Initial samples | 413 | 406 |
| Samples after QC | 379 (91.8%) | 377 (92.9%) |
| Initial SNPs | 36,901 | 700,000 |
| SNPs after mind/MAF/geno | 26,474 | 191,364 |
| SNPs after LD pruning | 1,187 (4.5%) | 12,975 (6.8%) |
| SNP retention from initial | 3.2% | 1.9% |
| Subpopulations represented | 6 | 6 |

### Table 2. Subpopulation Composition Before and After QC (44K)
| Subpopulation | 413 set | 379 set | Excluded | % change |
|---------------|---------|---------|----------|----------|
| ADMIX | 62 | 54 | 8 | −12.9% |
| AROMATIC | 14 | 13 | 1 | −7.1% |
| AUS | 57 | 49 | 8 | −14.0% |
| IND | 87 | 77 | 10 | −11.5% |
| TEJ | 96 | 93 | 3 | −3.1% |
| TRJ | 97 | 93 | 4 | −4.1% |

### Table 3. 44K Sensitivity Scenarios
| Dimension | Scenarios | Baseline |
|-----------|-----------|----------|
| Sample (S) | none, mind 0.10 (after QC), mind 0.10, mind 0.05, mind 0.02 | mind 0.10 |
| MAF (M) | none, 0.01, 0.03, 0.05, 0.10 | 0.05 |
| GENO (G) | none, 0.20, 0.10, 0.05, 0.02 | 0.05 |
| LD (L) | none, r² 0.8, r² 0.5, 50/5/0.2, 100/10/0.2 | 50/5/0.2 |

### Table 4. HDRA Sensitivity Scenarios
| Dimension | Scenarios | Baseline |
|-----------|-----------|----------|
| Sample (S) | none, mind 0.10*, mind 0.05 | mind 0.10 |
| MAF (M) | none, 0.01, 0.03, 0.05*, 0.10 | MAF 0.05 |
| GENO (G) | none, 0.20, 0.10, 0.05*, 0.02 | geno 0.05 |
| LD (L) | none, r² 0.8, r² 0.5, r² 0.2*, 100/10/0.2 | 50/5/0.2 |
\* Baseline scenario. S_002 (mind 0.02) was degenerate on HDRA and excluded.

### Table 5. FST Across Scenarios by Platform
| Scenario | 44K FST | HDRA FST |
|----------|---------|----------|
| Baseline | 0.423 | 0.326 |
| MAF none | 0.398 | 0.076 |
| MAF 0.01 | — | 0.120 |
| MAF 0.03 | — | 0.235 |
| MAF 0.10 | 0.443 | 0.443 |
| GENO none | 0.410 | 0.142 |
| GENO 0.02 | — | 0.350 |
| LD none | 0.609 | 0.512 |
| LD r² 0.5 | 0.434 | 0.354 |
| LD 100/10/0.2 | — | 0.299 |
| Sample none | 0.422 | 0.308 |
| Sample mind 0.05 | 0.291 | 0.324 |

### Table 6. ADMIXTURE Cross-Validation Errors (44K Baseline)
| K | Log-likelihood | CV error | ΔCV |
|---|---------------|----------|------|
| 1 | −464,095.94 | 1.04899 | — |
| 2 | −355,349.90 | 0.81137 | −0.23762 |
| 3 | −302,014.23 | 0.69623 | −0.11514 |
| 4 | −275,733.07 | 0.64413 | −0.05210 |
| 5 | −261,385.11 | 0.62135 | −0.02278 |
| 6 | −249,109.04 | 0.59605 | −0.02530 |
| 7 | −240,035.47 | 0.57871 | −0.01734 |
| 8 | −233,372.66 | 0.56911 | −0.00960 |
| 9 | −228,535.07 | 0.56226 | −0.00685 |
| 10 | −223,944.07 | 0.55569 | −0.00657 |
| 11 | −218,404.66 | 0.54761 | −0.00808 |
| 12 | −214,777.85 | 0.54633 | −0.00128 |

### Table 7. Cross-Platform Comparison Summary
| Metric | 44K | HDRA |
|--------|-----|------|
| Baseline samples | 379 | 377 |
| Baseline SNPs | 1,187 | 12,975 |
| Mean cluster agreement | 93.1% | 85.3% |
| Scenarios with ≥90% agreement | 81.0% | 64.7% |
| ARI range | 0.80–1.00 | 0.78–1.00 |
| Baseline FST | 0.423 | 0.326 |
| FST range (MAF dimension) | 0.398–0.443 | 0.076–0.443 |
| LD inflation factor | 1.44× | 1.81× |
| PC1 correlation range | 0.99–1.00 | 0.97–1.00 |
| He/PIC range | 0.27–0.31 | 0.059–0.370 |

### Table 8. DAPC Cluster Concordance with ADMIXTURE (44K Baseline, K = 5)
| | DAPC 1 | DAPC 2 | DAPC 3 | DAPC 4 | DAPC 5 | Total |
|---|:---:|:---:|:---:|:---:|:---:|:---:|
| ADMIXTURE cluster 1 | 0 | 87 | 0 | 0 | 0 | 87 |
| ADMIXTURE cluster 2 | 14 | 1 | 0 | 0 | 0 | 15 |
| ADMIXTURE cluster 3 | 0 | 1 | 50 | 0 | 0 | 51 |
| ADMIXTURE cluster 4 | 0 | 0 | 0 | 110 | 1 | 111 |
| ADMIXTURE cluster 5 | 0 | 0 | 0 | 0 | 115 | 115 |
| Total | 14 | 89 | 50 | 110 | 116 | 379 |

### Table 9. Kinship Matrix Summary (44K, 379 Accessions)
| Statistic | Value |
|-----------|-------|
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
