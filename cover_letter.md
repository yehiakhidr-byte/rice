[Date]

Dear Editor,

I am submitting our manuscript entitled "QC Parameter Sensitivity in Population-Genetic Inference: A Cross-Platform Comparison of the 44K and 700K SNP Arrays in the Rice Diversity Panel 1" for consideration as a research article in BMC Genomics.

**Methodological contribution.** The primary contribution of this work is a systematic, reproducible framework for evaluating how quality-control (QC) filtering decisions interact with genotyping platform density to shape population-genetic inferences. Using 34 QC scenarios across two array platforms (44K and 700K HDRA) in the well-characterized Rice Diversity Panel 1 (RDP1), we demonstrate that the relative importance of individual QC dimensions is platform-dependent: MAF filtering, which is the least sensitive dimension on the 44K array, becomes the dominant determinant of FST, genetic diversity, and cluster stability on the 700K HDRA platform (FST range: 0.076–0.443, a 5.8-fold difference). This finding has immediate practical implications as the field moves toward higher-density genotyping, and the stepwise sensitivity framework we present is directly transferable to other crop diversity panels.

**Reproducibility.** All data, code, and complete parameter documentation are publicly available in our GitHub repository (https://github.com/yehiakhidr-byte/rice). Specifically:

- All R and shell scripts for the complete analysis pipeline
- PLINK binary files for all 34 QC scenarios
- ADMIXTURE output (K = 1–10) for all scenarios
- Full supplementary tables including excluded-accession documentation (Table S14)
- Reproducibility documentation covering genotype file version, filtering history, software parameters, and bias assessment (REPRODUCIBILITY_DOCUMENTATION.md)

**Fit with BMC Genomics.** BMC Genomics' scope — "novel methods and techniques" and "comparative and evolutionary genomics" — aligns well with our work. The manuscript does not claim to discover new population structure in RDP1; rather, it asks a methodological question about how genotyping density affects the robustness of population-genetic metrics to QC variation. We believe this is the type of scientifically valid, methodologically rigorous contribution that BMC Genomics values.

Thank you for considering our manuscript. We look forward to your response.

Sincerely,
Yehia Khidr
[Affiliation]
[Email]
