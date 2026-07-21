# Rice Genomic Sensitivity Analysis

Complete reproducibility package for analysis of genetic diversity,
population structure, and QC sensitivity in the Rice Diversity Panel 1 (RDP1).

## Contents
- **R scripts**: Full analysis pipeline (QC, ADMIXTURE, PCA, FST, AMOVA, LD decay)
- **Shell scripts**: PLINK/ADMIXTURE batch processing
- **Documentation**: Complete parameter tables, excluded accession list, QC history
- **Metadata**: RDP1 passport data for all 413 accessions

## Reproducibility
See REPRODUCIBILITY_DOCUMENTATION.md for complete parameter documentation,
excluded accession list with metadata, and QC filtering history.

## Key results
- 16 sensitivity scenarios across 4 QC dimensions
- Cluster stability: 93.1% mean agreement after label alignment
- FST range: 0.29–0.61 (stable at 0.39–0.45 across reasonable QC choices)
- LD decay: r² ≈ 0.1 at ~388 kb (consistent with RDP1)
