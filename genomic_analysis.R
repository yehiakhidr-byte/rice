library(adegenet)
library(vcfR)
library(poppr)
library(RColorBrewer)
library(graphics)

setwd("/mnt/d/project genomic analysis")

# 1. ADMIXTURE ancestry bar plot (K=10)
q_mat <- as.matrix(read.table("rice_pruned_wsl.10.Q"))
fam <- read.table("rice_pruned_wsl.fam")
samples <- fam[, 2]
rownames(q_mat) <- samples

max_pop <- apply(q_mat, 1, which.max)
ord <- order(max_pop)
q_mat_ord <- q_mat[ord, ]

pdf("admixture_K10_barplot.pdf", width=12, height=4)
par(mar=c(4,4,1,1))
barplot(t(q_mat_ord), col=rainbow(10), border=NA,
        xlab="Individuals", ylab="Ancestry",
        main="ADMIXTURE Ancestry Proportions (K=10)")
dev.off()
cat("ADMIXTURE bar plot saved\n")

# CV error plot
cv_errors <- c(1.04161, 0.80398, 0.69746, 0.67526,
               0.61484, 0.60596, 0.59572, 0.57859,
               0.56743, 0.56397)
pdf("admixture_cv_errors.pdf", width=6, height=4)
plot(1:10, cv_errors, type="b", pch=19, col="blue",
     xlab="K", ylab="CV error",
     main="ADMIXTURE Cross-Validation Error")
dev.off()
cat("CV error plot saved\n")

# 2. DAPC
ped_file <- "/mnt/d/Rice accessions/RiceDiversity.44K.MSU6.Genotypes_PLINK/RiceDiversity_44K_Genotypes_PLINK/rice_pruned.raw"

if (file.exists(ped_file)) {
  obj <- read.PLINK(ped_file, parallel=FALSE)
  cat("Genind:", nInd(obj), "inds,", nLoc(obj), "loci\n")

  gt <- tab(obj)
  for (i in 1:ncol(gt)) {
    gt[is.na(gt[,i]), i] <- mean(gt[,i], na.rm=TRUE)
  }

  maxK <- min(20, nInd(obj) - 1)
  N <- nInd(obj)

  pca_res <- prcomp(gt, scale=FALSE, center=TRUE)
  n_pcs <- min(40, N - 1)
  pca_x <- pca_res$x[, 1:n_pcs]

  set.seed(42)
  grp <- find.clusters(pca_x, max.n.clust=maxK, n.iter=1e6, n.start=10,
                       stat="BIC", n.pca=n_pcs, n.clust=5)

  n_groups <- length(unique(grp$grp))
  cat("Clusters:", n_groups, "\n")

  n_da <- n_groups - 1
  dapc_npca <- min(30, N - 1, ncol(pca_x))
  dapc_res <- dapc(pca_x, grp$grp, n.pca=dapc_npca, n.da=n_da)

  pdf("dapc_scatter.pdf", width=7, height=7)
  scatter(dapc_res, scree.da=FALSE, legend=TRUE,
          posi.leg="bottomright", clab=0, cex=1.5)
  dev.off()
  cat("DAPC scatter plot saved\n")

  pdf("dapc_membership.pdf", width=10, height=4)
  matplot(dapc_res$posterior, type="l", lty=1, lwd=1,
          col=rainbow(n_groups),
          xlab="Individuals", ylab="Membership",
          main="DAPC Membership Probabilities")
  dev.off()
  cat("DAPC membership plot saved\n")
} else {
  cat("PLINK raw file not found.\n")
}

# 3. Kinship heatmap
kinship_file <- "/mnt/d/project genomic analysis/rice_kinship.rel"
if (file.exists(kinship_file)) {
  kins_raw <- read.table(kinship_file, fill=TRUE, col.names=1:413)
  kins <- as.matrix(kins_raw)
  cat("Kinship matrix:", nrow(kins), "x", ncol(kins), "\n")

  pdf("kinship_heatmap.pdf", width=10, height=9)
  kins[upper.tri(kins)] <- t(kins)[upper.tri(kins)]
  heatmap(kins, col=colorRampPalette(c("white","yellow","red"))(100),
          main="Kinship Matrix Heatmap", cexRow=0.5, cexCol=0.5)
  dev.off()
  cat("Kinship heatmap saved\n")

  pdf("kinship_distribution.pdf", width=6, height=4)
  hist(kins[upper.tri(kins)], breaks=50, col="steelblue",
       main="Pairwise Kinship Coefficients", xlab="Kinship coefficient")
  dev.off()
  cat("Kinship distribution saved\n")
}

cat("\nAll analyses complete!\n")