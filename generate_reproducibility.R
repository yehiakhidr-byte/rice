# ============================================================
# Reproducibility Documentation for Reviewer Response
# Generates: excluded accession table, QC history, parameter doc
# ============================================================
library(adegenet)

workdir <- "/home/yehia/sensitivity"
setwd(workdir)

# ---- 0. Load data ----
crossref <- read.csv("/mnt/d/project genomic analysis/RDP1_full_crossref.csv", 
                     stringsAsFactors=FALSE)
passport <- read.csv("/mnt/d/project genomic analysis/RDP1_413_complete_accession_passport.csv",
                     stringsAsFactors=FALSE)

fam_413 <- read.table("rice_pruned.fam", header=FALSE, stringsAsFactors=FALSE)
fam_409 <- read.table("rice_pruned2.fam", header=FALSE, stringsAsFactors=FALSE)
fam_379 <- read.table("rice_pruned3.fam", header=FALSE, stringsAsFactors=FALSE)

ids_413 <- fam_413$V2
ids_409 <- fam_409$V2
ids_379 <- fam_379$V2

# ---- 1. Excluded accessions ----
excl_1 <- setdiff(ids_413, ids_409)  # 4 removed by mind
excl_2 <- setdiff(ids_409, ids_379)  # 30 removed by stepwise pipeline
excl_all <- setdiff(ids_413, ids_379)

# Match numeric FAM IDs to NSFTV_ID
fam_to_nsftv <- function(fam_id) {
  paste0("NSFTV_", fam_id)
}

# Build exclusion table
build_exclusion_table <- function(fam_ids, step_label) {
  tbl <- data.frame(
    NSFTV_ID = fam_to_nsftv(fam_ids),
    FAM_ID = fam_ids,
    Exclusion_Step = step_label,
    stringsAsFactors=FALSE
  )
  # Merge with crossref
  merged <- merge(tbl, crossref, by="NSFTV_ID", all.x=TRUE, sort=FALSE)
  # Merge with passport for Subpopulation (may differ)
  merged2 <- merge(merged, passport[, c("NSFTV_ID", "Subpopulation")], 
                   by="NSFTV_ID", all.x=TRUE, sort=FALSE)
  names(merged2)[names(merged2) == "Subpopulation.x"] <- "Subpopulation_CrossRef"
  names(merged2)[names(merged2) == "Subpopulation.y"] <- "Subpopulation_Passport"
  merged2
}

cat("========================================\n")
cat("EXCLUDED ACCESSIONS: STEP 1 (mind, 413->409)\n")
cat("========================================\n")
t1 <- build_exclusion_table(excl_1, "Mind_413to409")
print(t1[, c("NSFTV_ID", "Accession_Name", "Country", "Subpopulation_Passport")])

cat("\n========================================\n")
cat("EXCLUDED ACCESSIONS: STEP 2 (pipeline, 409->379)\n")
cat("========================================\n")
t2 <- build_exclusion_table(excl_2, "Pipeline_409to379")
print(t2[, c("NSFTV_ID", "Accession_Name", "Country", "Subpopulation_Passport")])

cat("\n========================================\n")
cat("ALL 34 EXCLUDED ACCESSIONS\n")
cat("========================================\n")
t_all <- rbind(
  build_exclusion_table(excl_1, "Mind_413to409"),
  build_exclusion_table(excl_2, "Pipeline_409to379")
)
print(t_all[, c("NSFTV_ID", "Accession_Name", "Country", "Subpopulation_Passport", "Exclusion_Step")])

# Write table
write.csv(t_all, "/mnt/d/project genomic analysis/Table_S14_excluded_accessions.csv", row.names=FALSE)
cat("Table_S14_excluded_accessions.csv written\n")

# ---- 2. Subpopulation representation before/after ----
cat("\n========================================\n")
cat("SUBPOPULATION REPRESENTATION: 413 vs 379\n")
cat("========================================\n")

# Get subpopulation for all 413 accessions
merge_all <- merge(data.frame(NSFTV_ID = fam_to_nsftv(ids_413), stringsAsFactors=FALSE),
                   passport, by="NSFTV_ID", all.x=TRUE)
merge_379 <- merge(data.frame(NSFTV_ID = fam_to_nsftv(ids_379), stringsAsFactors=FALSE),
                   passport, by="NSFTV_ID", all.x=TRUE)

cat("413 set:\n")
print(table(merge_all$Subpopulation))
cat("\n379 set:\n")
print(table(merge_379$Subpopulation))
cat("\nDifference (excluded):\n")
diff_tbl <- table(merge_all$Subpopulation) - table(factor(merge_379$Subpopulation, levels=names(table(merge_all$Subpopulation))))
print(diff_tbl)

# ---- 3. Country representation change ----
tab_413_country <- sort(table(merge_all$Country), decreasing=TRUE)
tab_379_country <- sort(table(merge_379$Country), decreasing=TRUE)
cat("\nTop 10 countries in 413:\n")
print(head(tab_413_country, 10))
cat("\nTop 10 countries in 379:\n")
print(head(tab_379_country, 10))

# ---- 4. Excluded subpopulation counts ----
cat("\nExcluded by subpopulation:\n")
excl_subp <- merge(data.frame(NSFTV_ID = fam_to_nsftv(excl_all), stringsAsFactors=FALSE),
                   passport, by="NSFTV_ID", all.x=TRUE)
print(table(excl_subp$Subpopulation))

# ---- 5. Why the 4 accessions were excluded at 413->409 step ----
# Run PLINK missing report
cat("\n========================================\n")
cat("MISSING DATA FOR EXCLUDED ACCESSIONS\n")
cat("========================================\n")
# The mind filter removes samples with >10% missing genotypes
# PLINK generates a log file: rice_pruned.log
logfile <- "rice_pruned.log"
if (file.exists(logfile)) {
  cat("PLINK log for rice_pruned:\n")
  system(paste("grep -i 'missing\\|excluded\\|removed\\|people'", logfile))
}

# For the 409->379 step (30 accessions), these were removed during the
# stepwise MAF -> GENO -> LD pipeline. Let's check the actual reason.
# They could be due to:
#   - MAF 0.05 filtering removing informative SNPs for some samples
#   - GENO 0.05 removing SNPs with >5% missing
#   - Sample missingness creep after removing SNPs
# Actually, the 30 accessions are removed because the stepwise pipeline
# filters SNPs first, then re-checks sample missingness.
# Let's verify: the 409->379 step happens after all SNP filtering.

# Check rice_pruned2.log for mind filtering details
cat("\nPossible causes (from intermediate logs):\n")
for (logf in c("sens_D1_S3.log", "sens_D2_M3.log", "sens_D3_G3.log", "sens_D4_LD3.log")) {
  if (file.exists(logf)) {
    cat("  Log:", logf, "- ")
    system(paste("grep 'people.*remaining\\|removed.*people'", logf, "| head -1"))
  }
}

cat("\n===== REPRODUCIBILITY DOCUMENTATION GENERATED =====\n")
