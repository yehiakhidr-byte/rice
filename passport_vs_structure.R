library(adegenet)

setwd("/mnt/d/project genomic analysis")

# Load passport
passport <- read.csv("/mnt/d/Rice accessions/RDP1_full_crossref.csv", stringsAsFactors=FALSE)
passport$NSFTV_num <- as.numeric(gsub("NSFTV_", "", passport[,1]))

# Load fam
fam <- read.table("rice_pruned3.fam", stringsAsFactors=FALSE)
annot <- merge(fam, passport[, c("NSFTV_num", "Subpopulation", "Country", "Accession_Name")],
               by.x="V2", by.y="NSFTV_num", all.x=TRUE, sort=FALSE)

# Load ADMIXTURE Q matrices
q2 <- as.matrix(read.table("rice_pruned3.2.Q"))
q5 <- as.matrix(read.table("rice_pruned3.5.Q"))
annot$Admix2 <- apply(q2, 1, which.max)
annot$Admix5 <- apply(q5, 1, which.max)

# ========================================================
# TABLE A: ADMIXTURE K=2 vs PASSPORT SUBPOPULATION
# ========================================================
cat("============================================================\n")
cat("TABLE A: ADMIXTURE K=2 vs Passport Subpopulation\n")
cat("============================================================\n")
cross2 <- table(ADMIXTURE_K2=annot$Admix2, Subpopulation=annot$Subpopulation)
print(cross2)
cat("\nRow percentages:\n")
prop2 <- prop.table(cross2, 1) * 100
print(round(prop2, 1))

# ========================================================
# TABLE B: ADMIXTURE K=5 vs PASSPORT SUBPOPULATION
# ========================================================
cat("\n============================================================\n")
cat("TABLE B: ADMIXTURE K=5 vs Passport Subpopulation\n")
cat("============================================================\n")
cross5 <- table(ADMIXTURE_K5=annot$Admix5, Subpopulation=annot$Subpopulation)
print(cross5)
cat("\nRow percentages:\n")
prop5 <- prop.table(cross5, 1) * 100
print(round(prop5, 1))

cat("\nColumn percentages:\n")
col5 <- prop.table(cross5, 2) * 100
print(round(col5, 1))

# ========================================================
# TABLE C: PASSPORT SUBPOPULATION SUMMARY
# ========================================================
cat("\n============================================================\n")
cat("TABLE C: Passport Subpopulation Summary\n")
cat("============================================================\n")
cat(sprintf("%-15s %5s %s\n", "Subpopulation", "Count", "Countries"))
cat(strrep("-", 70), "\n")
for (sub in names(sort(table(annot$Subpopulation), decreasing=TRUE))) {
  sub_data <- annot[annot$Subpopulation == sub, ]
  countries <- names(sort(table(sub_data$Country), decreasing=TRUE))
  country_str <- paste(countries[1:min(5, length(countries))], collapse=", ")
  cat(sprintf("%-15s %5d  %s\n", sub, nrow(sub_data), country_str))
}

# ========================================================
# TABLE D: ADMIXTURE K=5 CLUSTER COMPOSITION
# ========================================================
cat("\n============================================================\n")
cat("TABLE D: ADMIXTURE K=5 Cluster Composition\n")
cat("============================================================\n")
cat(sprintf("%-10s %5s  %s\n", "Cluster", "Size", "Dominant subpopulation(s)"))
cat(strrep("-", 70), "\n")
for (cl in 1:5) {
  cl_data <- annot[annot$Admix5 == cl, ]
  sub_count <- sort(table(cl_data$Subpopulation), decreasing=TRUE)
  subs <- paste(names(sub_count), " (", sub_count, ")", sep="", collapse=", ")
  cat(sprintf("%-10d %5d  %s\n", cl, nrow(cl_data), subs))
}

# ========================================================
# TABLE E: TOP COUNTRIES PER SUBPOPULATION
# ========================================================
cat("\n============================================================\n")
cat("TABLE E: Top 5 Countries per Subpopulation\n")
cat("============================================================\n")
for (sub in c("IND", "AUS", "TEJ", "TRJ", "AROMATIC", "ADMIX")) {
  sub_data <- annot[annot$Subpopulation == sub, ]
  countries <- sort(table(sub_data$Country), decreasing=TRUE)[1:5]
  cat(sprintf("\n%s (n=%d):\n", sub, nrow(sub_data)))
  for (i in seq_along(countries)) {
    cat(sprintf("  %s (%d)\n", names(countries)[i], countries[i]))
  }
}

# ========================================================
# TABLE F: COUNTRY-LEVEL SUMMARY
# ========================================================
cat("\n============================================================\n")
cat("TABLE F: Country-level Subpopulation Distribution (top 10 countries)\n")
cat("============================================================\n")
country_table <- table(annot$Country, annot$Subpopulation)
top_ctry <- names(sort(rowSums(country_table), decreasing=TRUE)[1:10])
print(country_table[top_ctry, ])

cat("\nAll comparisons complete.\n")
