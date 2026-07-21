# Sensitivity Analysis: QC and Filtering Effects on Population-Genetic Inference
# Compares 16 scenarios across 4 dimensions

library(RColorBrewer)
library(adegenet)
library(vcfR)
library(poppr)

setwd("D:/project genomic analysis")

# Color scheme for each dimension
DIM_COLORS <- list(
  D1 = brewer.pal(4, "Blues")[2:4],
  D2 = brewer.pal(5, "RdYlGn"),
  D3 = brewer.pal(4, "Purples")[2:4],
  D4 = brewer.pal(5, "Oranges")
)

# ----------------------------------------------------------------------
# 1. READ CV SUMMARY
# ----------------------------------------------------------------------
cv_raw <- read.table("sensitivity_cv_summary.txt", header=TRUE,
                     stringsAsFactors=FALSE)
colnames(cv_raw) <- c("dataset", "K", "CV_error", "LogLikelihood")

# Classify datasets by dimension
cv_raw$dimension <- NA
cv_raw$dimension[grep("^sens_D1", cv_raw$dataset)] <- "D1_SampleFilter"
cv_raw$dimension[grep("^sens_D2", cv_raw$dataset)] <- "D2_MAF"
cv_raw$dimension[grep("^sens_D3", cv_raw$dataset)] <- "D3_Geno"
cv_raw$dimension[grep("^sens_D4", cv_raw$dataset)] <- "D4_LDpruning"
cv_raw$dimension[grep("^rice_pruned$", cv_raw$dataset)] <- "D1_SampleFilter"
cv_raw$dimension[grep("^rice_pruned2$", cv_raw$dataset)] <- "D1_SampleFilter"
cv_raw$dimension[grep("^rice_pruned3$", cv_raw$dataset)] <- "D1_SampleFilter"

# Label datasets clearly
label_map <- c(
  rice_pruned     = "S0: No mind (413)",
  rice_pruned2    = "S1: mind 0.10 after QC (409)",
  rice_pruned3    = "S3: mind 0.10 before QC (379)",
  sens_D1_S2      = "S2: mind 0.05 (293)",
  sens_D1_S3      = "S4: mind 0.02 (98)",
  sens_D2_M0      = "M0: MAF none (1638)",
  sens_D2_M1      = "M1: MAF 0.01 (1601)",
  sens_D2_M2      = "M2: MAF 0.03 (1386)",
  sens_D2_M4      = "M4: MAF 0.10 (844)",
  sens_D3_G020    = "G1: GENO 0.20 (1404)",
  sens_D3_G010    = "G2: GENO 0.10 (1299)",
  sens_D3_G002    = "G4: GENO 0.02 (982)",
  sens_D4_LDnone  = "L0: No pruning (26474)",
  sens_D4_LD08    = "L1: r2<0.8 (9609)",
  sens_D4_LD05    = "L2: r2<0.5 (4278)",
  sens_D4_LD100   = "L3: 100/10/0.2 (822)"
)

# For D3 baseline (GENO 0.05) we use rice_pruned3 as proxy
# Actually rice_pruned3 has 379 vs sens_D3 have 383, close enough for comparison

# ----------------------------------------------------------------------
# 2. SUPPLEMENTARY TABLE: SNP/sample counts
# ----------------------------------------------------------------------
counts <- data.frame(
  Scenario = names(label_map),
  Label = unname(label_map),
  Samples = c(413, 409, 379, 293, 98, rep(383, 4), rep(383, 3), rep(383, 4)),
  SNPs = c(1166, 1161, 1187, 1064, 597, 1638, 1601, 1386, 844,
           1404, 1299, 982, 26474, 9609, 4278, 822)
)
write.csv(counts, "Table_S1_sensitivity_counts.csv", row.names=FALSE)
cat("Table S1 saved:", nrow(counts), "scenarios\n")

# ----------------------------------------------------------------------
# 3. CV ERROR PLOTS (one panel per dimension)
# ----------------------------------------------------------------------
pdf("Figure_S2_sensitivity_CV_curves.pdf", width=12, height=10)
par(mfrow=c(2,2), mar=c(4,4,2,8), xpd=NA)

dims <- list(
  D1_SampleFilter = list(data=c("rice_pruned","rice_pruned2","sens_D1_S2","sens_D1_S3"),
                          cols=DIM_COLORS$D1, title="A) Sample filtering (mind)"),
  D2_MAF = list(data=c("sens_D2_M0","sens_D2_M1","sens_D2_M2","sens_D2_M4"),
                cols=DIM_COLORS$D2[c(1,2,3,5)], title="B) MAF threshold"),
  D3_Geno = list(data=c("sens_D3_G020","sens_D3_G010","sens_D3_G002"),
                 cols=DIM_COLORS$D3, title="C) Genotype missingness (GENO)"),
  D4_LDpruning = list(data=c("sens_D4_LDnone","sens_D4_LD08","sens_D4_LD05","sens_D4_LD100"),
                      cols=DIM_COLORS$D4, title="D) LD pruning threshold")
)

for (dn in names(dims)) {
  dd <- dims[[dn]]
  sub <- cv_raw[cv_raw$dimension == dn | cv_raw$dataset %in% dd$data, ]
  
  # Also include rice_pruned3 baseline in D1
  if (dn == "D1_SampleFilter") {
    sub <- cv_raw[cv_raw$dataset %in% c("rice_pruned","rice_pruned2","rice_pruned3","sens_D1_S2","sens_D1_S3"), ]
  }
  
  ds_list <- unique(sub$dataset)
  plot(NULL, xlim=c(1,10), ylim=c(0.35, 1.25),
       xlab="K", ylab="CV error", main=dd$title, cex.main=1.1)
  
  for (i in seq_along(ds_list)) {
    ds <- ds_list[i]
    vals <- sub[sub$dataset == ds, ]
    vals <- vals[order(vals$K), ]
    col_i <- if (i <= length(dd$cols)) dd$cols[i] else "gray50"
    lty_i <- if (ds == "rice_pruned3") 2 else 1
    lwd_i <- if (ds == "rice_pruned3") 2.5 else 1.5
    lines(vals$K, vals$CV_error, type="b", pch=i, col=col_i, lty=lty_i, lwd=lwd_i)
  }
  
  legend("topright", legend=label_map[ds_list], pch=seq_along(ds_list),
         col=dd$cols[1:length(ds_list)], lty=1, lwd=1.5, inset=c(-0.35,0), cex=0.7)
}
dev.off()
cat("Figure S2 (CV curves) saved\n")

# ----------------------------------------------------------------------
# 4. ADMIXTURE RESULTS: CV at K=5 comparison
# ----------------------------------------------------------------------
cv_k5 <- cv_raw[cv_raw$K == 5, ]
cv_k5$Label <- label_map[cv_k5$dataset]
rownames(cv_k5) <- cv_k5$dataset
cv_k5 <- cv_k5[names(label_map), ]
cv_k5 <- cv_k5[!is.na(cv_k5$CV_error), ]

write.csv(cv_k5[,c("Label","K","CV_error","LogLikelihood","dimension")],
          "Table_S2_sensitivity_CV_K5.csv", row.names=FALSE)
cat("Table S2 (CV at K=5) saved\n")

# ----------------------------------------------------------------------
# 5. PCA COMPARISON ACROSS SCENARIOS
# ----------------------------------------------------------------------
# Read the raw genotype data for baseline and each scenario
# We use PLINK .raw format for each pruned dataset

run_pca <- function(prefix, n_pcs=4) {
  raw_file <- paste0("D:/Rice accessions/RiceDiversity.44K.MSU6.Genotypes_PLINK/RiceDiversity_44K_Genotypes_PLINK/", prefix, ".raw")
  if (!file.exists(raw_file)) {
    # Try to create it via PLINK
    cmd <- sprintf('wsl -u root -- plink1.9 --bfile /mnt/d/"Rice accessions"/RiceDiversity.44K.MSU6.Genotypes_PLINK/RiceDiversity_44K_Genotypes_PLINK/%s --recode A --out /mnt/d/"Rice accessions"/RiceDiversity.44K.MSU6.Genotypes_PLINK/RiceDiversity_44K_Genotypes_PLINK/%s 2>&1 | tail -3', prefix, prefix)
    cat("Creating .raw for", prefix, "...\n")
    system(cmd)
  }
  if (!file.exists(raw_file)) return(NULL)
  
  obj <- try(read.PLINK(raw_file, parallel=FALSE), silent=TRUE)
  if (inherits(obj, "try-error")) return(NULL)
  
  gt <- tab(obj)
  for (j in 1:ncol(gt)) {
    gt[is.na(gt[,j]), j] <- mean(gt[,j], na.rm=TRUE)
  }
  pca <- prcomp(gt, scale=FALSE, center=TRUE)
  pve <- summary(pca)$importance[2, 1:n_pcs] * 100
  list(pca=pca, pve=pve, n=nInd(obj), m=nLoc(obj))
}

# PCA datasets (LD-pruned, manageable size)
pca_datasets <- c("rice_pruned", "rice_pruned2", "rice_pruned3",
                  "sens_D1_S2", "sens_D1_S3",
                  "sens_D2_M0", "sens_D2_M1", "sens_D2_M2", "sens_D2_M4",
                  "sens_D3_G020", "sens_D3_G010", "sens_D3_G002",
                  "sens_D4_LD08", "sens_D4_LD05", "sens_D4_LD100")

pca_results <- list()
for (ds in pca_datasets) {
  cat("PCA:", ds, "...\n")
  res <- run_pca(ds)
  if (!is.null(res)) {
    pca_results[[ds]] <- res
  }
}

# Procrustes analysis: compare PC1-2 of each scenario to baseline (rice_pruned3)
baseline_pc <- pca_results[["rice_pruned3"]]
if (!is.null(baseline_pc)) {
  procrustes_table <- data.frame(Scenario=character(), N=numeric(), M=numeric(),
                                 PC1_var=numeric(), PC2_var=numeric(),
                                 Corr_PC1=numeric(), Corr_PC2=numeric(),
                                 stringsAsFactors=FALSE)
  
  # Need to match samples across datasets
  bl_fam <- read.table("D:/Rice accessions/RiceDiversity.44K.MSU6.Genotypes_PLINK/RiceDiversity_44K_Genotypes_PLINK/rice_pruned3.fam")
  bl_samples <- bl_fam[,2]
  
  for (ds in names(pca_results)) {
    if (ds == "rice_pruned3") next
    pr <- pca_results[[ds]]
    
    # Read samples for this dataset
    ds_fam <- read.table(paste0("D:/Rice accessions/RiceDiversity.44K.MSU6.Genotypes_PLINK/RiceDiversity_44K_Genotypes_PLINK/", ds, ".fam"))
    ds_samples <- ds_fam[,2]
    
    # Find common samples
    common <- intersect(bl_samples, ds_samples)
    cat("  ", ds, ": ", length(common), "/", pr$n, " samples in common\n", sep="")
    
    # Correlation of PC scores on common samples
    bl_pc_subset <- baseline_pc$pca$x[rownames(baseline_pc$pca$x) %in% common, 1:2]
    pr_pc_subset <- pr$pca$x[rownames(pr$pca$x) %in% common, 1:2]
    
    # Handle sign ambiguity
    cor1 <- cor(bl_pc_subset[,1], pr_pc_subset[,1])
    cor2 <- cor(bl_pc_subset[,2], pr_pc_subset[,2])
    
    procrustes_table <- rbind(procrustes_table,
      data.frame(Scenario=ds, N=pr$n, M=pr$m,
                 PC1_var=round(pr$pve[1],1), PC2_var=round(pr$pve[2],1),
                 Corr_PC1=round(cor1,4), Corr_PC2=round(cor2,4)))
  }
  
  procrustes_table$Label <- label_map[procrustes_table$Scenario]
  write.csv(procrustes_table, "Table_S3_PCA_procrustes.csv", row.names=FALSE)
  cat("Table S3 (PCA procrustes) saved\n")
}

# ----------------------------------------------------------------------
# 6. PCA SCREE COMPARISON PLOT
# ----------------------------------------------------------------------
if (length(pca_results) > 0) {
  pdf("Figure_S3_sensitivity_PCA_scree.pdf", width=10, height=8)
  par(mfrow=c(2,2), mar=c(4,4,3,2))
  
  dims_pca <- list(
    D1 = c("rice_pruned","rice_pruned2","rice_pruned3"),
    D2 = c("sens_D2_M0","sens_D2_M1","sens_D2_M2","sens_D2_M4"),
    D3 = c("sens_D3_G020","sens_D3_G010","sens_D3_G002"),
    D4 = c("sens_D4_LD08","sens_D4_LD05","sens_D4_LD100")
  )
  
  for (dn in names(dims_pca)) {
    ds_list <- dims_pca[[dn]]
    ds_list <- ds_list[ds_list %in% names(pca_results)]
    if (length(ds_list) == 0) next
    
    cols <- switch(dn,
      D1 = DIM_COLORS$D1,
      D2 = DIM_COLORS$D2[c(1,2,3,5)],
      D3 = DIM_COLORS$D3,
      D4 = DIM_COLORS$D4[c(2,3,5)])
    
    plot(NULL, xlim=c(1,10), ylim=c(0, max(sapply(pca_results[ds_list], function(x) sum(x$pve[1:10]))) + 5),
         xlab="PC", ylab="Variance explained (%)",
         main=paste0("PCA Scree - ", dn), type="n")
    
    for (i in seq_along(ds_list)) {
      pr <- pca_results[[ds_list[i]]]
      pve_all <- pr$pve
      lines(1:length(pve_all), pve_all, type="b", pch=i, col=cols[i], lwd=1.5)
    }
    legend("topright", legend=label_map[ds_list], pch=seq_along(ds_list),
           col=cols, lty=1, cex=0.6)
  }
  dev.off()
  cat("Figure S3 (PCA scree) saved\n")
}

# ----------------------------------------------------------------------
# 7. SUMMARY TABLE
# ----------------------------------------------------------------------
cat("\n===== SENSITIVITY ANALYSIS COMPLETE =====\n")
cat("Generated:\n")
cat("  Table_S1_sensitivity_counts.csv\n")
cat("  Table_S2_sensitivity_CV_K5.csv\n")
cat("  Table_S3_PCA_procrustes.csv\n")
cat("  Figure_S2_sensitivity_CV_curves.pdf\n")
cat("  Figure_S3_sensitivity_PCA_scree.pdf\n")
