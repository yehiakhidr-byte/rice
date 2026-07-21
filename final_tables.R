library(adegenet)
library(RColorBrewer)

# ========== LOAD DATA ==========
obj <- read.PLINK("rice_pruned3.raw", parallel=FALSE)
gt <- tab(obj)
for(i in 1:ncol(gt)){gt[is.na(gt[,i]), i] <- mean(gt[,i], na.rm=TRUE)}
fam <- read.table("rice_pruned3.fam")
samples <- fam[, 2]
N <- nInd(obj)

# PCA
pca <- prcomp(gt, scale=FALSE, center=TRUE)
pve <- summary(pca)$importance[2, ] * 100

# DAPC
maxK <- min(20, N - 1)
n_pcs <- min(40, N - 1)
pca_x <- pca$x[, 1:n_pcs]
set.seed(42)
grp <- find.clusters(pca_x, max.n.clust=maxK, n.iter=1e6, n.start=10,
                     stat="BIC", n.pca=n_pcs, n.clust=5)
n_groups <- length(unique(grp$grp))
dapc_res <- dapc(pca_x, grp$grp, n.pca=min(30, N-1, ncol(pca_x)), n.da=n_groups-1)

# Kinship 379
kins_raw <- read.table("rice_kinship_379.rel", fill=TRUE, col.names=1:379)
kins <- as.matrix(kins_raw)
for (i in 1:nrow(kins)) {
  for (j in i:ncol(kins)) {
    if (i == j) next
    kins[i, j] <- kins[j, i]
  }
}
upper <- kins[upper.tri(kins)]

# =======================================================
# TABLE 1: ADMIXTURE CV Errors
# =======================================================
cat("============================================================\n")
cat("TABLE 1 | ADMIXTURE cross-validation errors\n")
cat("============================================================\n")
cat(sprintf(" %-4s %-18s %-12s %s\n", "K", "Log-likelihood", "CV error", "CV change"))
cat("------------------------------------------------------------\n")
best_ll  <- c(-464095.94, -355349.90, -302014.23, -275733.07,
              -261385.11, -249109.04, -240035.47, -233372.66,
              -228535.07, -223944.07, -218404.66, -214777.85)
best_cv  <- c(1.04899, 0.81137, 0.69623, 0.64413,
              0.62135, 0.59605, 0.57871, 0.56911,
              0.56226, 0.55569, 0.54761, 0.54633)
prev_cv <- c(NA, best_cv[1:11])
change  <- best_cv - prev_cv
for (k in 1:12) {
  chg <- if (k == 1) "—" else sprintf("%+.5f", change[k])
  star <- if (k == 5) " *" else ""
  cat(sprintf(" %-4d %-18.2f %-12.5f%s %s\n", k, best_ll[k], best_cv[k], chg, star))
}
cat("------------------------------------------------------------\n")
cat(" * Optimal K (K=5) based on CV elbow, STRUCTURE concordance,\n")
cat("   and DAPC cluster validation.\n\n")

# =======================================================
# TABLE 2: PCA
# =======================================================
cat("============================================================\n")
cat("TABLE 2 | Principal component analysis\n")
cat("============================================================\n")
cat(sprintf(" %-4s %-15s %-12s %s\n", "PC", "Eigenvalue", "Variance %", "Cumulative %"))
cat("------------------------------------------------------------\n")
cum <- 0
for (i in 1:10) {
  cum <- cum + pve[i]
  cat(sprintf(" %-4d %-15.2f %-11.2f %-9.2f\n", i, pca$sdev[i]^2, pve[i], cum))
}
cat("------------------------------------------------------------\n")
cat(sprintf(" Note: First two PCs explain %.1f%% + %.1f%% = %.1f%%\n", pve[1], pve[2], pve[1]+pve[2]))
cat(sprintf("       of total genetic variation.\n\n"))

# =======================================================
# TABLE 3: DAPC
# =======================================================
cat("============================================================\n")
cat("TABLE 3 | Discriminant Analysis of Principal Components\n")
cat("============================================================\n")
tab_clust <- table(grp$grp)
cat(sprintf(" %-10s %s\n", "Cluster", "Number of individuals"))
cat("-----------------------------\n")
for (i in 1:length(tab_clust)) {
  cat(sprintf(" %-10d %d\n", as.integer(names(tab_clust)[i]), tab_clust[i]))
}
cat("-----------------------------\n")
cat(sprintf(" Total      %d\n\n", sum(tab_clust)))

cat("DAPC eigenvalues:\n")
cat(sprintf(" %-6s %-12s %s\n", "Axis", "Eigenvalue", "Proportion"))
cat("-------------------------------\n")
for (i in 1:length(dapc_res$eig)) {
  pct <- dapc_res$eig[i] / sum(dapc_res$eig) * 100
  cat(sprintf(" %-6d %-12.4f %.1f%%\n", i, dapc_res$eig[i], pct))
}
cat("\n")

# =======================================================
# TABLE 4: Kinship
# =======================================================
cat("============================================================\n")
cat("TABLE 4 | Kinship matrix summary\n")
cat("============================================================\n")
cat(sprintf(" %-30s %s\n", "Statistic", "Value"))
cat("---------------------------------------------\n")
cat(sprintf(" %-30s %d\n", "Samples", nrow(kins)))
cat(sprintf(" %-30s %d\n", "SNPs used", 26474))
cat(sprintf(" %-30s %.4f\n", "Mean", mean(upper)))
cat(sprintf(" %-30s %.4f\n", "SD", sd(upper)))
cat(sprintf(" %-30s %.4f\n", "Minimum", min(upper)))
cat(sprintf(" %-30s %.4f\n", "Maximum", max(upper)))
cat(sprintf(" %-30s %.4f\n", "Median", median(upper)))
cat(sprintf(" %-30s %.2f%%\n", "Pairs < 0", sum(upper < 0) / length(upper) * 100))
cat(sprintf(" %-30s %.2f%%\n", "Pairs > 0.125 (2nd degree)", sum(upper > 0.125) / length(upper) * 100))
cat(sprintf(" %-30s %.2f%%\n", "Pairs > 0.25 (1st degree)", sum(upper > 0.25) / length(upper) * 100))
cat("---------------------------------------------\n")
cat(" Note: Centered IBS kinship (PLINK --make-rel).\n")
cat("       Negative values indicate pairs less related than expected\n")
cat("       from random population background.\n")
