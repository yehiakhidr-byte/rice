library(adegenet)

# Read data
obj <- read.PLINK("rice_pruned3.raw", parallel=FALSE)
gt <- tab(obj)
for(i in 1:ncol(gt)){gt[is.na(gt[,i]), i] <- mean(gt[,i], na.rm=TRUE)}

fam <- read.table("rice_pruned3.fam")
samples <- fam[, 2]

# PCA
pca <- prcomp(gt, scale=FALSE, center=TRUE)
pve <- summary(pca)$importance[2, ] * 100

# DAPC (same params as genomic_analysis4.R)
N <- nInd(obj)
maxK <- min(20, N - 1)
n_pcs <- min(40, N - 1)
pca_x <- pca$x[, 1:n_pcs]

set.seed(42)
grp <- find.clusters(pca_x, max.n.clust=maxK, n.iter=1e6, n.start=10,
                     stat="BIC", n.pca=n_pcs, n.clust=5)
n_groups <- length(unique(grp$grp))
n_da <- n_groups - 1
dapc_npca <- min(30, N - 1, ncol(pca_x))
dapc_res <- dapc(pca_x, grp$grp, n.pca=dapc_npca, n.da=n_da)

# === TABLE 1: ADMIXTURE CV errors ===
cat("=== Table 1: ADMIXTURE Cross-Validation Errors ===\n")
cat(sprintf("%-6s %-15s %-15s\n", "K", "Best LL", "CV Error"))
cat(strrep("-", 40), "\n")
best_cv <- c(1.04899, 0.81137, 0.69623, 0.64413,
             0.62135, 0.59605, 0.57871, 0.56911,
             0.56226, 0.55569, 0.54761, 0.54633)
best_ll <- c(-464095.94, -355349.90, -302014.23, -275733.07,
             -261385.11, -249109.04, -240035.47, -233372.66,
             -228535.07, -223944.07, -218404.66, -214777.85)
for (k in 1:12) {
  flag <- if (k == 5) " <-- optimal" else ""
  cat(sprintf("%-6d %-15.2f %-15.5f%s\n", k, best_ll[k], best_cv[k], flag))
}

# === TABLE 2: PCA variance explained (top 10 PCs) ===
cat("\n=== Table 2: PCA Variance Explained (top 10 PCs) ===\n")
cat(sprintf("%-6s %-15s %-15s\n", "PC", "Eigenvalue", "Variance (%)"))
cat(strrep("-", 40), "\n")
for (i in 1:10) {
  cat(sprintf("%-6d %-15.2f %-15.2f\n", i, pca$sdev[i]^2, pve[i]))
}

# === TABLE 3: DAPC cluster sizes ===
cat("\n=== Table 3: DAPC Cluster Sizes ===\n")
tab_clust <- table(grp$grp)
cat(sprintf("%-15s %-10s\n", "Cluster", "N"))
cat(strrep("-", 30), "\n")
for (i in 1:length(tab_clust)) {
  cat(sprintf("%-15d %-10d\n", as.integer(names(tab_clust)[i]), tab_clust[i]))
}

cat(sprintf("\nTotal: %d individuals\n", sum(tab_clust)))

# DAPC eigenvalues
cat("\n=== Table 4: DAPC Eigenvalues ===\n")
cat(sprintf("%-15s %-15s\n", "Axis", "Eigenvalue"))
cat(strrep("-", 30), "\n")
for (i in 1:length(dapc_res$eig)) {
  pct <- dapc_res$eig[i] / sum(dapc_res$eig) * 100
  cat(sprintf("%-15d %-15.4f (%.1f%%)\n", i, dapc_res$eig[i], pct))
}

# === TABLE 5: Kinship summary ===
cat("\n=== Table 5: Kinship Summary ===\n")
kins_raw <- read.table("rice_kinship.rel", fill=TRUE, col.names=1:413)
kins <- as.matrix(kins_raw)
kins[upper.tri(kins)] <- t(kins)[upper.tri(kins)]
upper <- kins[upper.tri(kins)]
cat(sprintf("%-25s %s\n", "Statistic", "Value"))
cat(strrep("-", 40), "\n")
cat(sprintf("%-25s %d x %d\n", "Matrix dimensions", nrow(kins), ncol(kins)))
cat(sprintf("%-25s %.4f\n", "Mean", mean(upper)))
cat(sprintf("%-25s %.4f\n", "SD", sd(upper)))
cat(sprintf("%-25s %.4f\n", "Min", min(upper)))
cat(sprintf("%-25s %.4f\n", "Max", max(upper)))
cat(sprintf("%-25s %.4f\n", "Median", median(upper)))
cat(sprintf("%-25s %.2f%%\n", "Negative values", sum(upper < 0) / length(upper) * 100))
cat(sprintf("%-25s %.2f%%\n", ">0.25 (close kin)", sum(upper > 0.25) / length(upper) * 100))
cat(sprintf("%-25s %.2f%%\n", ">0.125 (2nd degree)", sum(upper > 0.125) / length(upper) * 100))

cat("\nAll tables generated.\n")
