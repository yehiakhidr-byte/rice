# HDRA RDP1 Sensitivity Analysis
# Comparison with 44K RDP1 results
# Run after hdra_run.sh pipeline completes

library(ggplot2)
library(reshape2)
library(dplyr)
library(vegan)    # for ARI/NMI
library(ape)      # for AMOVA

setwd("/home/yehia/hdra_rdp1")

# Scenario tags (matching the bash pipeline)
scenarios <- c(
  "hdra_S_none", "hdra_S_010", "hdra_S_005", "hdra_S_002",
  "hdra_M_none", "hdra_M_001", "hdra_M_003", "hdra_M_005", "hdra_M_010",
  "hdra_G_none", "hdra_G_020", "hdra_G_010", "hdra_G_005", "hdra_G_002",
  "hdra_L_none", "hdra_L_08", "hdra_L_05", "hdra_L_02", "hdra_L_100"
)

short_names <- c(
  "S_none","S_mind10","S_mind5","S_mind2",
  "M_none","M_maf1","M_maf3","M_maf5","M_maf10",
  "G_none","G_geno20","G_geno10","G_geno5","G_geno2",
  "L_none","L_ld08","L_ld05","L_ld02","L_ld100"
)

dimensions <- c(rep("Sample\nfiltering",4), rep("MAF\nthreshold",5),
                rep("GENO\nthreshold",5), rep("LD\npruning",5))

# Load sample data
fam <- read.table("hdra_rdp1_raw.fam", stringsAsFactors=FALSE)
n_samples <- nrow(fam)

# Load subpopulation map
subpop <- read.table("subpop_map.txt", stringsAsFactors=FALSE)
colnames(subpop) <- c("hash", "subpop")
pop_map <- setNames(subpop$subpop, subpop$hash)

# =========================================================
# 1. Read scenario stats and ADMIXTURE results
# =========================================================
stats <- read.csv("s.csv", stringsAsFactors=FALSE)
cat("Scenario stats:\n")
print(stats)

# Read K=5 Q-files for all scenarios
q5_list <- list()
for (i in seq_along(scenarios)) {
  s <- scenarios[i]
  fn <- paste0(s, "_K5.Q")
  if (file.exists(fn)) {
    q5_list[[short_names[i]]] <- as.matrix(read.table(fn))
  }
}

# Assign clusters based on max Q
cluster_assignments <- sapply(q5_list, function(q) apply(q, 1, which.max))
colnames(cluster_assignments) <- names(q5_list)

# =========================================================
# 2. Cluster label alignment (greedy)
# =========================================================
align_clusters_greedy <- function(target, query) {
  k <- max(target, query)
  mapping <- integer(k)
  used <- logical(k)
  for (ti in 1:k) {
    best_val <- -1
    best_qi <- 0
    for (qi in 1:k) {
      if (!used[qi]) {
        agreement <- sum(target == ti & query == qi)
        if (agreement > best_val) {
          best_val <- agreement
          best_qi <- qi
        }
      }
    }
    mapping[best_qi] <- ti
    used[best_qi] <- TRUE
  }
  return(mapping)
}

# Baseline: L_02 (LD r2=0.2, corresponds to baseline scenario)
baseline <- "L_02"
baseline_clust <- cluster_assignments[, baseline]

aligned <- matrix(NA, nrow=n_samples, ncol=length(scenarios))
colnames(aligned) <- short_names
aligned[, baseline] <- baseline_clust

for (s in setdiff(short_names, baseline)) {
  q <- cluster_assignments[, s]
  k <- max(q, baseline_clust)
  map <- align_clusters_greedy(baseline_clust, q)
  aligned[, s] <- map[q]
}

# =========================================================
# 3. Cluster stability
# =========================================================
agreement <- sapply(setdiff(short_names, baseline), function(s) {
  sum(aligned[, s] == aligned[, baseline]) / n_samples
})
cat("\nCluster agreement with baseline (L_02):\n")
print(round(agreement, 4))

cat("\nMean agreement:", round(mean(agreement), 4))
cat("\nProportion with >=0.9 agreement:", round(mean(agreement >= 0.9), 4))

# =========================================================
# 4. Confusion matrices
# =========================================================
pdf("hdra_confusion_matrices.pdf", width=12, height=8)
par(mfrow=c(4,5), mar=c(3,3,3,1))
for (s in short_names) {
  tab <- table(aligned[, baseline], aligned[, s])
  image(t(tab), col=heat.colors(12), axes=FALSE, main=s, cex.main=0.7)
  axis(1, at=seq(0,1,length.out=nrow(tab))-0.5/nrow(tab), labels=1:nrow(tab), cex.axis=0.6)
  axis(2, at=seq(0,1,length.out=ncol(tab))-0.5/ncol(tab), labels=1:ncol(tab), cex.axis=0.6)
}
dev.off()
cat("Confusion matrices saved.\n")

# =========================================================
# 5. Pairwise ARI / NMI
# =========================================================
calc_ari <- function(c1, c2) {
  library(mclust)
  adjustedRandIndex(c1, c2)
}
calc_nmi <- function(c1, c2) {
  library(aricode)
  NMI(c1, c2)
}

ari_mat <- matrix(NA, length(short_names), length(short_names))
rownames(ari_mat) <- colnames(ari_mat) <- short_names
nmi_mat <- ari_mat

for (i in 1:length(short_names)) {
  for (j in 1:length(short_names)) {
    ari_mat[i,j] <- calc_ari(aligned[,i], aligned[,j])
    nmi_mat[i,j] <- calc_nmi(aligned[,i], aligned[,j])
  }
}

write.csv(ari_mat, "hdra_ari_matrix.csv")
write.csv(nmi_mat, "hdra_nmi_matrix.csv")

pdf("hdra_ari_nmi_heatmaps.pdf", width=10, height=8)
par(mfrow=c(1,2), mar=c(8,8,3,3))
image(1:ncol(ari_mat), 1:nrow(ari_mat), t(ari_mat[nrow(ari_mat):1,]), 
      col=heat.colors(20), axes=FALSE, xlab="", ylab="", main="ARI")
axis(1, at=1:length(short_names), labels=short_names, las=2, cex.axis=0.6)
axis(2, at=length(short_names):1, labels=short_names, las=2, cex.axis=0.6)
image(1:ncol(nmi_mat), 1:nrow(nmi_mat), t(nmi_mat[nrow(nmi_mat):1,]), 
      col=heat.colors(20), axes=FALSE, xlab="", ylab="", main="NMI")
axis(1, at=1:length(short_names), labels=short_names, las=2, cex.axis=0.6)
axis(2, at=length(short_names):1, labels=short_names, las=2, cex.axis=0.6)
dev.off()

# =========================================================
# 6. PCA correlations
# =========================================================
pca_list <- list()
for (s in scenarios) {
  fn <- paste0(s, "_pca.eigenvec")
  if (file.exists(fn)) {
    pca_list[[s]] <- read.table(fn, header=FALSE)[,3:12]
  }
}

# Correlate PCs between scenarios
pca_cor_list <- list()
for (i in seq_along(scenarios)) {
  s1 <- scenarios[i]
  if (!s1 %in% names(pca_list)) next
  cors <- numeric(length(scenarios))
  for (j in seq_along(scenarios)) {
    s2 <- scenarios[j]
    if (!s2 %in% names(pca_list)) next
    # PC1 vs PC1
    cors[j] <- abs(cor(pca_list[[s1]][,1], pca_list[[s2]][,1]))
  }
  names(cors) <- short_names
  pca_cor_list[[short_names[i]]] <- cors
}
pca_cor_mat <- do.call(rbind, pca_cor_list)
write.csv(pca_cor_mat, "hdra_pca_correlations.csv")

pdf("hdra_pca_cor_heatmap.pdf", width=8, height=7)
image(1:ncol(pca_cor_mat), 1:nrow(pca_cor_mat), t(pca_cor_mat[nrow(pca_cor_mat):1,]),
      col=heat.colors(20), axes=FALSE, main="PC1 Correlation across scenarios")
axis(1, at=1:ncol(pca_cor_mat), labels=colnames(pca_cor_mat), las=2, cex.axis=0.6)
axis(2, at=nrow(pca_cor_mat):1, labels=rownames(pca_cor_mat), las=2, cex.axis=0.6)
dev.off()

# =========================================================
# 7. Genetic diversity
# =========================================================
div_list <- list()
for (s in scenarios) {
  fn <- paste0(s, "_hwe.hwe")
  if (!file.exists(fn)) next
  hwe <- try(read.table(fn, header=TRUE, skip=1))
  if (inherits(hwe, "try-error")) next
  # He = expected heterozygosity
  he <- mean(2 * hwe$O.HET. * (1 - hwe$O.HET.), na.rm=TRUE)
  ho <- mean(hwe$O.HET., na.rm=TRUE)
  pic <- mean(1 - (hwe$O.HET.^2 + (1-hwe$O.HET.)^2), na.rm=TRUE)
  div_list[[s]] <- c(He=he, Ho=ho, PIC=pic)
}
div_df <- do.call(rbind, div_list)
rownames(div_df) <- short_names[1:nrow(div_df)]
write.csv(div_df, "hdra_diversity.csv")
cat("\nDiversity metrics:\n")
print(round(div_df, 4))

# Diversity plot
pdf("hdra_diversity.pdf", width=10, height=5)
par(mfrow=c(1,3), mar=c(8,4,3,1))
for (m in c("He", "Ho", "PIC")) {
  barplot(div_df[,m], names.arg=rownames(div_df), las=2, cex.names=0.6,
          main=m, ylim=c(0, max(div_df[,m], na.rm=TRUE)*1.2))
}
dev.off()

# =========================================================
# 8. FST comparison
# =========================================================
fst_list <- list()
for (s in scenarios) {
  fn <- paste0(s, "_fst.fst")
  if (file.exists(fn)) {
    fst <- try(read.table(fn, header=TRUE))
    if (!inherits(fst, "try-error") && nrow(fst) > 0) {
      fst_list[[s]] <- mean(fst[,ncol(fst)], na.rm=TRUE)
    }
  }
}
if (length(fst_list) > 0) {
  fst_vec <- unlist(fst_list)
  write.csv(data.frame(scenario=names(fst_vec), FST=fst_vec), "hdra_fst.csv", row.names=FALSE)
  cat("\nFST values:\n")
  print(round(fst_vec, 4))
}

# =========================================================
# 9. ADMIXTURE CV error for K selection
# =========================================================
cv_list <- list()
for (s in scenarios) {
  for (K in 1:10) {
    fn <- paste0(s, "_K", K, ".log")
    if (file.exists(fn)) {
      log_content <- readLines(fn)
      cv_line <- grep("CV", log_content, value=TRUE)
      if (length(cv_line) > 0) {
        cv_val <- as.numeric(gsub(".*CV error \\(cross-validation\\):? ?", "", cv_line))
        cv_list[[length(cv_list)+1]] <- data.frame(scenario=s, K=K, CV=cv_val)
      }
    }
  }
}
if (length(cv_list) > 0) {
  cv_df <- do.call(rbind, cv_list)
  write.csv(cv_df, "hdra_cv_errors.csv", row.names=FALSE)
}

# =========================================================
# 10. Summary table
# =========================================================
summary_df <- data.frame(
  Scenario = short_names,
  Dimension = dimensions,
  Samples = stats$samples[match(scenarios, stats$scenario)],
  SNPs = stats$snps[match(scenarios, stats$scenario)],
  Agreement = NA,
  Mean_He = div_df[short_names, "He"],
  Mean_Ho = div_df[short_names, "Ho"],
  Mean_PIC = div_df[short_names, "PIC"]
)
summary_df$Agreement[match(names(agreement), summary_df$Scenario)] <- agreement

write.csv(summary_df, "hdra_summary.csv", row.names=FALSE)
cat("\n=== HDRA Analysis Complete ===\n")
print(summary_df[,1:5])

# =========================================================
# 11. Compare with 44K results
# =========================================================
cat("\n\n=== Comparison with 44K Analysis ===\n")

# Load 44K summary
if (file.exists("/home/yehia/sensitivity/summary.csv")) {
  k44 <- read.csv("/home/yehia/sensitivity/summary.csv", stringsAsFactors=FALSE)
  cat("44K samples:", unique(k44$Samples)[1], "| HDRA RDP1 samples:", n_samples, "\n")
  cat("44K baseline SNPs:", subset(k44, Scenario=="S3_baseline")$SNPs, "| HDRA baseline SNPs:", stats$snps[stats$scenario==scenarios[which(short_names=="L_02")]], "\n")
  cat("44K mean agreement:", round(mean(k44$Agreement, na.rm=TRUE), 4), "| HDRA:", round(mean(agreement, na.rm=TRUE), 4), "\n")
}

cat("\n\nDone.\n")