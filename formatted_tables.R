# Generate formatted Tables 5-8 as text and PDF

# ---- TABLE 5: ADMIXTURE CV Errors ----
table5_text <- "
Table 5. ADMIXTURE cross-validation errors for K = 1-12.
+-----+------------------+-----------+-----------+
|  K  | Log-likelihood   | CV error  | Delta CV  |
+-----+------------------+-----------+-----------+
|  1  | -464,095.94      |   1.04899 |     —     |
|  2  | -355,349.90      |   0.81137 |  -0.23762 |
|  3  | -302,014.23      |   0.69623 |  -0.11514 |
|  4  | -275,733.07      |   0.64413 |  -0.05210 |
|  5* | -261,385.11      |   0.62135 |  -0.02278 |
|  6  | -249,109.04      |   0.59605 |  -0.02530 |
|  7  | -240,035.47      |   0.57871 |  -0.01734 |
|  8  | -233,372.66      |   0.56911 |  -0.00960 |
|  9  | -228,535.07      |   0.56226 |  -0.00685 |
| 10  | -223,944.07      |   0.55569 |  -0.00657 |
| 11  | -218,404.66      |   0.54761 |  -0.00808 |
| 12  | -214,777.85      |   0.54633 |  -0.00128 |
+-----+------------------+-----------+-----------+
Footnotes: Values are from the best replicate (highest log-likelihood)
out of 10 independent runs per K. K = 5* is the optimal K based on
CV error elbow, concordance with STRUCTURE, and DAPC validation.
"
cat(table5_text)
writeLines(table5_text, "Table5_admixture_cv.txt")

# ---- TABLE 6: Concordance ----
table6_text <- "
Table 6. Concordance between ADMIXTURE and DAPC cluster assignments.
(A) K = 2

               DAPC cluster
ADMIXTURE      Cluster 1    Cluster 2    Total
-----------------------------------------------
Cluster 1            6          140       146
Cluster 2          233            0       233
-----------------------------------------------
Total              239          140       379

(B) K = 5

               DAPC cluster
ADMIXTURE      1     2     3     4     5   Total
--------------------------------------------------
Cluster 1      0    87     0     0     0      87
Cluster 2     14     1     0     0     0      15
Cluster 3      0     1    50     0     0      51
Cluster 4      0     0     0   110     1     111
Cluster 5      0     0     0     0   115     115
--------------------------------------------------
Total         14    89    50   110   116     379

Footnotes: Values are counts of individuals assigned to each
combination of methods. At K = 5, clusters are numbered by
ADMIXTURE assignment order; corresponding DAPC clusters may
have different numeric labels.
"
cat(table6_text)
writeLines(table6_text, "Table6_concordance.txt")

# ---- TABLE 7: DAPC ----
table7_text <- "
Table 7. DAPC cluster sizes and eigenvalues.
(A) Cluster sizes

               Number of
Cluster        individuals
--------------------------
    1              14
    2              89
    3              50
    4             110
    5             116
--------------------------
  Total           379

(B) Discriminant analysis eigenvalues

         Eigen-     Variance
Axis     value      explained (%)
----------------------------------
  1    7,284.61        61.8
  2    2,478.59        21.0
  3    1,068.33         9.1
  4      955.24         8.1
----------------------------------
  Total               100.0
"
cat(table7_text)
writeLines(table7_text, "Table7_dapc.txt")

# ---- TABLE 8: Kinship ----
table8_text <- "
Table 8. Summary statistics of the centered IBS kinship matrix
for 379 rice accessions (26,474 QC-passing SNPs).

Statistic                              Value
----------------------------------------------
Samples                                   379
SNPs used                              26,474
Mean                                   -0.0054
Standard deviation                      0.8087
Minimum                                -1.0512
Maximum                                 3.6250
Median                                  0.0670
Proportion of pairs < 0                48.17%
Proportion of pairs > 0.125            48.65%
  (second-degree relatives)
Proportion of pairs > 0.25             39.66%
  (first-degree relatives)
----------------------------------------------
Footnotes: Centered IBS kinship calculated
using PLINK --make-rel with the full set of
26,474 QC-passing markers. Negative values
indicate pairs less related than expected
under random population background.
"
cat(table8_text)
writeLines(table8_text, "Table8_kinship.txt")

cat("\nAll formatted tables saved as text files.\n")
