library(adegenet)
library(vcfR)
library(poppr)
library(RColorBrewer)
library(graphics)

setwd("/mnt/d/project genomic analysis")

# 1. ADMIXTURE CV error plot (best replicate per K)
best_cv <- c(1.04899, 0.81137, 0.69623, 0.64413,
             0.62135, 0.59605, 0.57871, 0.56911,
             0.56226, 0.55569, 0.54761, 0.54633)
pdf("admixture_cv_errors_v4.pdf", width=6, height=4)
plot(1:12, best_cv, type="b", pch=19, col="blue",
     xlab="K", ylab="CV error",
     main="ADMIXTURE CV Error (best of 10 reps, 379 samples)")
points(5, best_cv[5], pch=19, col="red", cex=1.5)
text(5, best_cv[5], "K=5", pos=3, col="red", font=2)
dev.off()
cat("CV error plot saved\n")

# 2a. ADMIXTURE main bar plot (K=5, matches STRUCTURE and DAPC)
q5 <- as.matrix(read.table("rice_pruned3.5.Q"))
fam <- read.table("rice_pruned3.fam")
samples <- fam[, 2]
rownames(q5) <- samples

max_pop5 <- apply(q5, 1, which.max)
ord5 <- order(max_pop5)
q5_ord <- q5[ord5, ]

pdf("admixture_main_barplot_v4.pdf", width=12, height=4)
par(mar=c(4,4,1,1))
barplot(t(q5_ord), col=rainbow(5), border=NA,
        xlab="Individuals", ylab="Ancestry",
        main="ADMIXTURE Ancestry Proportions (K=5, 379 samples)")
dev.off()
cat("ADMIXTURE main bar plot (K=5) saved\n")

# 2b. ADMIXTURE supplementary bar plot (K=12)
q12 <- as.matrix(read.table("rice_pruned3.12.Q"))
rownames(q12) <- samples

max_pop12 <- apply(q12, 1, which.max)
ord12 <- order(max_pop12)
q12_ord <- q12[ord12, ]

pdf("admixture_supp_barplot_v4.pdf", width=12, height=4)
par(mar=c(4,4,1,1))
barplot(t(q12_ord), col=rainbow(12), border=NA,
        xlab="Individuals", ylab="Ancestry",
        main="ADMIXTURE Ancestry Proportions (K=12, Supplementary)")
dev.off()
cat("ADMIXTURE supp bar plot (K=12) saved\n")

# 3. PCA from PLINK binary data
ped_file <- "/mnt/d/project genomic analysis/rice_pruned3.raw"
if (file.exists(ped_file)) {
  obj <- read.PLINK(ped_file, parallel=FALSE)
  cat("Genind:", nInd(obj), "inds,", nLoc(obj), "loci\n")

  gt <- tab(obj)
  for (i in 1:ncol(gt)) {
    gt[is.na(gt[,i]), i] <- mean(gt[,i], na.rm=TRUE)
  }

  pca_res <- prcomp(gt, scale=FALSE, center=TRUE)
  pve <- summary(pca_res)$importance[2, 1:4] * 100

  pdf("pca_scatter_v4.pdf", width=7, height=6)
  plot(pca_res$x[,1], pca_res$x[,2], pch=19, cex=0.6,
       col="steelblue",
       xlab=paste0("PC1 (", round(pve[1], 1), "%)"),
       ylab=paste0("PC2 (", round(pve[2], 1), "%)"),
       main="PCA of 379 rice samples (1,187 LD-pruned SNPs)")
  dev.off()
  cat("PCA scatter plot saved\n")

  pdf("pca_pve_v4.pdf", width=6, height=4)
  plot(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, type="b", pch=19,
       col="darkgreen", xlab="PC", ylab="Variance explained (%)",
       main="PCA Scree Plot")
  dev.off()
  cat("PCA scree plot saved\n")

  # 4. DAPC
  N <- nInd(obj)
  maxK <- min(20, N - 1)
  n_pcs <- min(40, N - 1)
  pca_x <- pca_res$x[, 1:n_pcs]

  set.seed(42)
  grp <- find.clusters(pca_x, max.n.clust=maxK, n.iter=1e6, n.start=10,
                       stat="BIC", n.pca=n_pcs, n.clust=5)

  n_groups <- length(unique(grp$grp))
  cat("DAPC clusters:", n_groups, "\n")

  n_da <- n_groups - 1
  dapc_npca <- min(30, N - 1, ncol(pca_x))
  dapc_res <- dapc(pca_x, grp$grp, n.pca=dapc_npca, n.da=n_da)

  pdf("dapc_scatter_v4.pdf", width=7, height=7)
  scatter(dapc_res, scree.da=FALSE, legend=TRUE,
          posi.leg="bottomright", clab=0, cex=1.5)
  dev.off()
  cat("DAPC scatter plot saved\n")

  pdf("dapc_membership_v4.pdf", width=10, height=4)
  matplot(dapc_res$posterior, type="l", lty=1, lwd=1,
          col=rainbow(n_groups),
          xlab="Individuals", ylab="Membership",
          main="DAPC Membership Probabilities")
  dev.off()
  cat("DAPC membership plot saved\n")
}

# 5. Kinship heatmap
kinship_file <- "/mnt/d/project genomic analysis/rice_kinship.rel"
if (file.exists(kinship_file)) {
  kins_raw <- read.table(kinship_file, fill=TRUE, col.names=1:413)
  kins <- as.matrix(kins_raw)
  cat("Kinship matrix:", nrow(kins), "x", ncol(kins), "\n")

  kins[upper.tri(kins)] <- t(kins)[upper.tri(kins)]
  pdf("kinship_heatmap_v4.pdf", width=10, height=9)
  heatmap(kins, col=colorRampPalette(c("white","yellow","red"))(100),
          main="Kinship Matrix Heatmap", cexRow=0.5, cexCol=0.5)
  dev.off()
  cat("Kinship heatmap saved\n")

  pdf("kinship_distribution_v4.pdf", width=6, height=4)
  hist(kins[upper.tri(kins)], breaks=50, col="steelblue",
       main="Pairwise Kinship Coefficients", xlab="Kinship coefficient")
  dev.off()
  cat("Kinship distribution saved\n")
}

cat("\nAll v4 analyses complete!\n")
