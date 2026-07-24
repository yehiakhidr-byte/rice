# HDRA RDP1 Sensitivity Analysis
setwd("/home/yehia/hdra_rdp1")

library(mclust)   # for adjustedRandIndex
library(aricode)  # for NMI

# Scenario definitions
scenarios <- c(
  "hdra_S_none","hdra_S_010","hdra_S_005",
  "hdra_M_none","hdra_M_001","hdra_M_003","hdra_M_005","hdra_M_010",
  "hdra_G_none","hdra_G_020","hdra_G_010","hdra_G_005","hdra_G_002",
  "hdra_L_none","hdra_L_08","hdra_L_05","hdra_L_02","hdra_L_100"
)

short_names <- c(
  "S_none","S_mind10","S_mind5",
  "M_none","M_maf1","M_maf3","M_maf5","M_maf10",
  "G_none","G_geno20","G_geno10","G_geno5","G_geno2",
  "L_none","L_ld08","L_ld05","L_ld02","L_ld100"
)

dimensions <- c(rep("Sample\nfiltering",3), rep("MAF\nthreshold",5),
                rep("GENO\nthreshold",5), rep("LD\npruning",5))

# Load subpopulation map
subpop <- read.table("subpop_map.txt", stringsAsFactors=FALSE)
colnames(subpop) <- c("hash", "subpop")

# Track which scenarios share which sample sets
fam_files <- sapply(scenarios, function(s) {
  if (file.exists(paste0(s, ".pruned.fam"))) paste0(s, ".pruned.fam")
  else if (file.exists(paste0(s, ".fam"))) paste0(s, ".fam")
  else NA
})

# Get sample IDs per scenario
sample_ids <- lapply(seq_along(scenarios), function(i) {
  if (is.na(fam_files[i])) return(NULL)
  read.table(fam_files[i], stringsAsFactors=FALSE)[,2]
})
names(sample_ids) <- short_names

# Baseline: L_02 (same as 44K baseline)
baseline_name <- "L_ld02"
baseline_idx <- which(short_names == baseline_name)
baseline_ids <- sample_ids[[baseline_idx]]

cat("=== HDRA Sensitivity Analysis ===\n")
cat("Samples per scenario:\n")
print(sapply(sample_ids, length))

# Intersection of all samples (except S_none which has more)
common_samples <- Reduce(intersect, sample_ids[sapply(sample_ids, length) == length(baseline_ids)])
cat("\nSamples in intersection (all except S_none):", length(common_samples), "\n")

# =========================================================
# 1. Read K=5 ADMIXTURE results
# =========================================================
q5 <- list()
for (i in seq_along(scenarios)) {
  s <- scenarios[i]
  fn <- paste0(s, "_K5.Q")
  if (file.exists(fn)) {
    q5[[short_names[i]]] <- as.matrix(read.table(fn))
  }
}

# 2. Greedy cluster alignment
align_clusters <- function(target, query) {
  k <- max(target, query)
  map <- integer(k)
  used <- logical(k)
  for (t in 1:k) {
    best_val <- -1
    best_q <- 0
    for (q in 1:k) {
      if (!used[q]) {
        agree <- sum(target == t & query == q)
        if (agree > best_val) {
          best_val <- agree
          best_q <- q
        }
      }
    }
    map[best_q] <- t
    used[best_q] <- TRUE
  }
  return(map)
}

# For scenarios with same sample count as baseline, align directly
aligned <- list()
baseline_q <- q5[[baseline_name]]
baseline_clust <- apply(baseline_q, 1, which.max)

for (nm in names(q5)) {
  q <- q5[[nm]]
  clust <- apply(q, 1, which.max)
  
  if (length(clust) == length(baseline_clust)) {
    map <- align_clusters(baseline_clust, clust)
    aligned[[nm]] <- map[clust]
  } else {
    # Different sample count - use common IDs
    ids <- sample_ids[[nm]]
    base_ids <- sample_ids[[baseline_name]]
    common <- intersect(ids, base_ids)
    idx_this <- match(common, ids)
    idx_base <- match(common, base_ids)
    
    if (length(common) > 10) {
      map <- align_clusters(baseline_clust[idx_base], clust[idx_this])
      aligned[[nm]] <- map[clust]
    } else {
      aligned[[nm]] <- clust
    }
  }
}

# 3. Cluster agreement with baseline
cat("\n=== Cluster Agreement with L_ld02 ===\n")
agreement <- c()
for (nm in setdiff(names(aligned), baseline_name)) {
  if (length(aligned[[nm]]) == length(aligned[[baseline_name]])) {
    agree <- sum(aligned[[nm]] == aligned[[baseline_name]]) / length(aligned[[nm]])
  } else {
    common <- intersect(sample_ids[[nm]], sample_ids[[baseline_name]])
    idx_nm <- match(common, sample_ids[[nm]])
    idx_bl <- match(common, sample_ids[[baseline_name]])
    agree <- sum(aligned[[nm]][idx_nm] == aligned[[baseline_name]][idx_bl]) / length(common)
  }
  agreement[nm] <- round(agree, 4)
}
print(agreement)

cat("\nMean agreement:", round(mean(agreement), 4), "\n")
cat("Agreement >= 0.9:", round(mean(agreement >= 0.9) * 100, 1), "%\n")

# 4. Confusion matrices
pdf("hdra_confusion_matrices.pdf", width=12, height=8)
par(mfrow=c(3,6), mar=c(3,3,3,1))
for (nm in short_names) {
  if (length(aligned[[baseline_name]]) == length(aligned[[nm]])) {
    tab <- table(aligned[[baseline_name]], aligned[[nm]])
    image(t(tab), col=terrain.colors(12), axes=FALSE, main=nm, cex.main=0.6)
    if (nrow(tab) > 1) axis(1, at=seq(0,1,l=nrow(tab)), labels=1:nrow(tab), cex.axis=0.5)
    if (ncol(tab) > 1) axis(2, at=seq(0,1,l=ncol(tab)), labels=1:ncol(tab), cex.axis=0.5)
  }
}
dev.off()
cat("Confusion matrices PDF saved.\n")

# 5. Pairwise ARI/NMI
cat("\n=== Pairwise ARI/NMI ===\n")
n <- length(short_names)
ari_mat <- matrix(NA, n, n)
nmi_mat <- matrix(NA, n, n)
rownames(ari_mat) <- colnames(ari_mat) <- short_names
rownames(nmi_mat) <- colnames(nmi_mat) <- short_names

for (i in 1:n) {
  for (j in 1:n) {
    if (length(aligned[[short_names[i]]]) == length(aligned[[short_names[j]]])) {
      ari_mat[i,j] <- adjustedRandIndex(aligned[[short_names[i]]], aligned[[short_names[j]]])
      nmi_mat[i,j] <- NMI(aligned[[short_names[i]]], aligned[[short_names[j]]])
    }
  }
}

write.csv(ari_mat, "hdra_ari.csv")
write.csv(nmi_mat, "hdra_nmi.csv")
cat("ARI range:", round(range(ari_mat, na.rm=TRUE), 3), "\n")
cat("NMI range:", round(range(nmi_mat, na.rm=TRUE), 3), "\n")

# 6. PCA correlations
cat("\n=== PCA Correlations ===\n")
pca_cor <- matrix(NA, length(scenarios), length(scenarios))
rownames(pca_cor) <- colnames(pca_cor) <- short_names

for (i in seq_along(scenarios)) {
  s1 <- scenarios[i]
  fn1 <- paste0(s1, "_pca.eigenvec")
  if (!file.exists(fn1)) next
  pc1 <- try(read.table(fn1, header=FALSE)[,3])
  if (inherits(pc1, "try-error")) next
  
  for (j in seq_along(scenarios)) {
    s2 <- scenarios[j]
    fn2 <- paste0(s2, "_pca.eigenvec")
    if (!file.exists(fn2)) next
    pc2 <- try(read.table(fn2, header=FALSE)[,3])
    if (inherits(pc2, "try-error")) next
    
    if (length(pc1) == length(pc2)) {
      pca_cor[i,j] <- abs(cor(pc1, pc2))
    }
  }
}

write.csv(pca_cor, "hdra_pca_cor.csv")
cat("PC1 correlation range:", round(range(pca_cor, na.rm=TRUE), 3), "\n")

# 7. Diversity metrics (from HWE)
cat("\n=== Diversity Metrics ===\n")
div_list <- list()
for (i in seq_along(scenarios)) {
  s <- scenarios[i]
  fn <- paste0(s, "_hwe.hwe")
  if (!file.exists(fn)) next
  hwe <- try(read.table(fn, header=TRUE, comment.char=""))
  if (inherits(hwe, "try-error") || nrow(hwe) < 10) next
  
  ho <- mean(hwe$O.HET., na.rm=TRUE)
  he <- mean(hwe$E.HET., na.rm=TRUE)
  pic <- mean(1 - ((1-hwe$O.HET.)^2 + hwe$O.HET.^2), na.rm=TRUE)
  
  div_list[[short_names[i]]] <- c(He=he, Ho=ho, PIC=pic)
}
div_df <- do.call(rbind, div_list)
write.csv(div_df, "hdra_diversity.csv")
print(round(div_df, 4))

# 8. CV errors
cat("\n=== CV Errors for K selection ===\n")
cv_list <- list()
for (s in scenarios) {
  for (K in 1:10) {
    logfn <- paste0(s, "_K", K, ".log")
    if (file.exists(logfn)) {
      # ADMIXTURE log file has CV error in it
      logtxt <- readLines(logfn)
      cv_line <- grep("CV error", logtxt, value=TRUE)
      if (length(cv_line) > 0) {
        cv_val <- as.numeric(gsub(".*CV error.*?:? ?", "", cv_line))
        cv_list[[length(cv_list)+1]] <- data.frame(scenario=s, K=K, CV=cv_val)
      }
    }
  }
}
cv_df <- do.call(rbind, cv_list)
write.csv(cv_df, "hdra_cv_errors.csv", row.names=FALSE)

# 9. Summary table
cat("\n=== Summary Table ===\n")
stats <- read.csv("s.csv", stringsAsFactors=FALSE)
summary_df <- data.frame(
  Scenario = short_names,
  Dimension = dimensions,
  Samples = stats$samples[match(scenarios, stats$scenario)],
  SNPs = stats$snps[match(scenarios, stats$scenario)],
  Agreement = agreement[short_names],
  He = div_df[short_names, "He"],
  Ho = div_df[short_names, "Ho"],
  PIC = div_df[short_names, "PIC"]
)
rownames(summary_df) <- NULL
write.csv(summary_df, "hdra_summary.csv", row.names=FALSE)
print(summary_df[,1:5])

# =========================================================
# 10. Compare with 44K
# =========================================================
cat("\n\n=== Comparison with 44K Results ===\n")
if (file.exists("/home/yehia/sensitivity/summary.csv")) {
  k44 <- read.csv("/home/yehia/sensitivity/summary.csv", stringsAsFactors=FALSE)
  cat("44K samples:", unique(k44$Samples)[1], "\n")
  cat("HDRA RDP1 samples (baseline):", stats$samples[which(short_names=="L_ld02")], "\n\n")
  
  cat("44K baseline SNPs:", subset(k44, Scenario=="S3_baseline")$SNPs, "\n")
  cat("HDRA baseline SNPs:", stats$snps[which(short_names=="L_ld02")], "\n\n")
  
  cat("44K mean cluster agreement:", round(mean(k44$Agreement, na.rm=TRUE)*100, 1), "%\n")
  cat("HDRA mean cluster agreement:", round(mean(agreement, na.rm=TRUE)*100, 1), "%\n\n")
  
  # FST comparison
  hdra_fst <- read.csv("hdra_fst_summary.csv", stringsAsFactors=FALSE)
  cat("44K FST range:", round(range(k44$FST, na.rm=TRUE), 3), "\n")
  cat("HDRA FST range:", round(range(hdra_fst$Mean_FST, na.rm=TRUE), 3), "\n")
  cat("44K baseline FST:", round(subset(k44, Scenario=="S3_baseline")$FST, 4), "\n")
  cat("HDRA baseline FST:", round(hdra_fst$Mean_FST[which(short_names=="L_ld02")], 4), "\n")
}

cat("\n=== HDRA Analysis Complete ===\n")