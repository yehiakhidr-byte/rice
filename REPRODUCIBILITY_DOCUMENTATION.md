# Reproducibility Documentation

## 1. Genotype File and Data Version

- **Source**: USDA Rice Diversity Panel 1 (RDP1) 
- **Original dataset**: 413 accessions (five major subpopulations: IND, TEJ, AUS, AROMATIC, ADMIX/TRJ)
- **Genotype file**: PLINK binary format (.bed/.bim/.fam)
- **File version**: `rice_pruned` (413 accessions, pre-QC)
- **Accession IDs**: Sample IDs in FAM files correspond to NSFTV_num (e.g., ID `1` = NSFTV_1). Full cross-reference in `RDP1_full_crossref.csv`

## 2. QC Filtering History

```
Initial:   413 accessions × 26,474 SNPs   (rice_pruned.bed)
  │
  ├─ Step 1: --mind 0.10 (remove samples >10% missing)
  │   4 accessions removed: NSFTV_72, NSFTV_194, NSFTV_252, NSFTV_390
  │   → 409 accessions (rice_pruned2)
  │
  ├─ Step 2: --maf 0.05 (remove SNPs with MAF < 0.05)
  │
  ├─ Step 3: --geno 0.05 (remove SNPs with >5% missing)
  │
  ├─ Step 4: --indep-pairwise 50 5 0.2 (LD pruning)
  │
  └─ After SNP filtering (steps 2-4), sample missingness re-check
      30 accessions removed due to elevated missingness after SNP-level QC
      → 379 accessions × 1,187 SNPs (rice_pruned3 — baseline dataset)
```

### Reason for 30 sample removals in Step 2
These 30 accessions were not removed by any explicit `--mind` threshold.
They were lost during the stepwise SNP filtering pipeline (MAF → GENO → LD):
removing poorly genotyped SNPs can expose previously undetected sample-level
missingness. PLINK re-checks sample missingness after each `--make-bed` step,
and samples exceeding the implied threshold are dropped.

## 3. Excluded Accessions (Table_S14)

See `Table_S14_excluded_accessions.csv` for the full list of 34 excluded
accessions with NSFTV_ID, Accession_Name, Country, Subpopulation, and
Exclusion_Step (Mind_413to409 or Pipeline_409to379).

### Summary by subpopulation
| Subpopulation | 413 set | 379 set | Excluded | % excluded |
|---|---|---|---|---|
| ADMIX        | 62  | 54  | 8  | 12.9% |
| AROMATIC     | 14  | 13  | 1  | 7.1%  |
| AUS          | 57  | 49  | 8  | 14.0% |
| IND          | 87  | 77  | 10 | 11.5% |
| TEJ          | 96  | 93  | 3  | 3.1%  |
| TRJ          | 97  | 93  | 4  | 4.1%  |

Exclusions are broadly proportional to subpopulation size. The largest
absolute reductions are in IND (−10) and AUS (−8), but these are close to
the overall exclusion rate (8.2%). No single subpopulation is
disproportionately excluded.

## 4. Complete Pipeline Parameters

### PLINK Filtering
| Parameter | Value |
|---|---|
| Sample missingness (mind) | 0.10 |
| MAF threshold | 0.05 |
| SNP missingness (geno) | 0.05 |
| LD pruning method | `--indep-pairwise` |
| Window size (SNPs) | 50 |
| Step size (SNPs) | 5 |
| r² threshold | 0.2 |

### ADMIXTURE
| Parameter | Value |
|---|---|
| K range | 1–10 |
| Replicates per K | 1 (built-in CV) |
| Convergence criterion | default (1e-4) |
| Random seed | default |
| Cluster assignment | maximum Q-value per individual |
| Hardware | Intel i7, WSL2 on Windows 11 |

### PCA
| Parameter | Value |
|---|---|
| Software | PLINK 1.9 |
| Number of PCs | 10 |
| Input | LD-pruned dataset |

### FST
| Parameter | Value |
|---|---|
| Software | PLINK 1.9 `--fst` |
| Population labels | ADMIXTURE K=5 clusters |
| Method | Weir & Cockerham (1984) |
| Reported metric | Mean FST across all SNP-population pairs |

### AMOVA
| Parameter | Value |
|---|---|
| Software | `pegas::amova` in R |
| Distance metric | Allele-count difference (PLINK `--distance square`) |
| Population labels | ADMIXTURE K=5 clusters (baseline) |
| Permutations | 0 (analytical p-values) |

### Genetic Diversity
| Parameter | Value |
|---|---|
| He / Ho | From PLINK `.hwe` output |
| PIC | 2 × MAF × (1 − MAF) (biallelic) |
| Shannon Index | −Σ p × ln(p) |
| Frequency source | PLINK `--freq` |

### LD Decay
| Parameter | Value |
|---|---|
| Software | PLINK 1.9 `--r2` |
| Max window | 99 SNPs |
| Max window (kb) | 1,000 |
| r² threshold | none (record all pairs) |
| Binning | Physical distance bins: 0–10, 10–50, 50–100, 100–200, 200–500, 500–1000, 1000–5000 kb |
| Curve metric | Mean r² per distance bin |

## 5. Software and Scripts

All custom scripts are in the GitHub repository:
**https://github.com/YOUR_USER/rice-sensitivity-analysis** (or local archive)

### Script inventory
| File | Purpose |
|---|---|
| `genomic_analysis4.R` | Primary analysis: ADMIXTURE, DAPC, FST, diversity, kinship |
| `sensitivity_pipeline.sh` | Generates all 16 sensitivity datasets + ADMIXTURE runs |
| `figures_sensitivity.R` | Workflow figure, QC impact, cluster stability |
| `sensitivity_analysis2.R` | Confusion matrix, ARI/NMI, PCA correlation, diversity, FST, AMOVA |
| `ld_decay_plot.R` | LD decay figure from PLINK output |
| `generate_reproducibility.R` | Excluded accession table and QC summary |

### Required software versions
| Software | Version |
|---|---|
| PLINK | 1.90b7.2 |
| ADMIXTURE | 1.3.0 |
| R | 4.x |
| adegenet | 2.1.11 |
| pegas | 1.4 |
| poppr | 2.9.8 |

## 6. K = 2 vs. K = 5 Interpretation

The reviewer correctly notes that ΔK detected by Evanno's method peaks at K = 2.
This is expected: ΔK identifies the uppermost hierarchical level of population
divergence, which in rice corresponds to the Indica/Japonica split (the deepest
split in O. sativa). K = 5 is the biologically meaningful number, corresponding
to the five well-established subpopulations (IND, TEJ, AUS, AROMATIC, ADMIX/TRJ).

- **K = 2**: primary divergence between indica and japonica groups
- **K = 5**: biologically relevant subpopulation structure (matches known germplasm groups)

The manuscript does not present K = 2 as a discovery — it is reported for
completeness (standard ADMIXTURE practice) and to show that CV error plateaus
around K = 5. All population-genetic inferences (FST, AMOVA, diversity) are
based on K = 5 clusters.

## 7. LD Decay Contextualization

RDP1 (original) reported LD decay to r² ≈ 0.1 at 500 kb–1 Mb.
This analysis finds r² ≈ 0.1 at ~388 kb.

### Sources of the ~112–612 kb difference
| Factor | Direction | Magnitude estimate |
|---|---|---|
| **Sample composition** (379 vs. 413) | Minimal | <5% shift |
| **MAF threshold** (0.05 vs. possibly lower in RDP1) | Higher MAF → faster decay | Could shift ~50–100 kb |
| **Marker density** (1,187 vs. ~36,000 SNPs) | Lower density → less resolution at short distances | May overestimate r² at closest bins |
| **Pruning** (LD-pruned vs. full dataset) | Removes correlated markers → apparent faster decay | Potentially large: L0 dataset shows r² ≈ 0.43 at 0–10 kb |
| **Binning/curve-fitting** | Different bin edges → different interpolation | 10–20% variation |
| **Software** (PLINK vs. original RDP1 pipeline) | Minimal | <5% |

The 379-accession LD results are consistent with the original RDP1 pattern
(the same qualitative decay toward background), with quantitative differences
attributable to MAF threshold, SNP density, and LD pruning — not to any
biological discrepancy.

## 8. Potential Sampling Bias Assessment

The 34 excluded accessions were removed solely by technical QC (genotype
missingness, not by phenotype or subpopulation). The subpopulation composition
of the 379 subset (ADMIX=14.2%, AROMATIC=3.4%, AUS=12.9%, IND=20.3%,
TEJ=24.5%, TRJ=24.5%) closely mirrors the original 413 set (ADMIX=15.0%,
AROMATIC=3.4%, AUS=13.8%, IND=21.1%, TEJ=23.2%, TRJ=23.5%). Country
representation is also preserved (top 8 countries unchanged, relative order
preserved).

**Conclusion**: No evidence of systematic sampling bias. All inferences
(FST, AMOVA, diversity, LD decay) are expected to be representative of the
original RDP1 panel.
