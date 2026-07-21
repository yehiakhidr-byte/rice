library(adegenet)
library(RColorBrewer)
library(grDevices)

setwd("/mnt/d/project genomic analysis")

# ========== LOAD DATA ==========
obj <- read.PLINK("rice_pruned3.raw", parallel=FALSE)
gt <- tab(obj)
for(i in 1:ncol(gt)){gt[is.na(gt[,i]), i] <- mean(gt[,i], na.rm=TRUE)}
fam <- read.table("rice_pruned3.fam")
samples <- fam[, 2]
N <- nInd(obj)

pca <- prcomp(gt, scale=FALSE, center=TRUE)
pve <- summary(pca)$importance[2, ] * 100
n_pcs <- min(40, N - 1)
pca_x <- pca$x[, 1:n_pcs]

# ADMIXTURE Q matrices
q2 <- as.matrix(read.table("rice_pruned3.2.Q"))
rownames(q2) <- samples
q5 <- as.matrix(read.table("rice_pruned3.5.Q"))
rownames(q5) <- samples

# DAPC K=2
set.seed(42)
grp2 <- find.clusters(pca_x, max.n.clust=20, n.iter=1e6, n.start=10,
                      stat="BIC", n.pca=n_pcs, n.clust=2)
dapc2 <- dapc(pca_x, grp2$grp, n.pca=min(30, N-1, ncol(pca_x)), n.da=1)

# DAPC K=5
set.seed(42)
grp5 <- find.clusters(pca_x, max.n.clust=20, n.iter=1e6, n.start=10,
                      stat="BIC", n.pca=n_pcs, n.clust=5)
dapc5 <- dapc(pca_x, grp5$grp, n.pca=min(30, N-1, ncol(pca_x)), n.da=4)

# ============================================================
# FIGURE 4: ADMIXTURE multi-panel
# CV error (a) + K=2 barplot (b) + K=5 barplot (c)
# ============================================================
best_cv <- c(1.04899, 0.81137, 0.69623, 0.64413,
             0.62135, 0.59605, 0.57871, 0.56911,
             0.56226, 0.55569, 0.54761, 0.54633)

pdf("Figure4_admixture.pdf", width=10, height=9)
layout(matrix(c(1,2,3), nrow=3))
par(mar=c(4,4,2,1))

# Panel (a): CV error
plot(1:12, best_cv, type="b", pch=19, col="blue", cex=1.2,
     xlab="K", ylab="CV error", main="(a) Cross-validation error",
     cex.main=1.2, cex.lab=1.1)
points(5, best_cv[5], pch=19, col="red", cex=2)
text(5, best_cv[5], "K=5", pos=3, col="red", font=2, cex=1.2)
box(lwd=1)

# Panel (b): K=2 barplot
par(mar=c(4,2,2,1))
max_pop2 <- apply(q2, 1, which.max)
ord2 <- order(max_pop2)
barplot(t(q2[ord2, ]), col=c("tomato","steelblue"), border=NA,
        xlab="Individuals", ylab="Ancestry", main="(b) K = 2",
        cex.main=1.2, cex.lab=1.1)

# Panel (c): K=5 barplot
par(mar=c(4,2,2,1))
max_pop5 <- apply(q5, 1, which.max)
ord5 <- order(max_pop5)
barplot(t(q5[ord5, ]), col=rainbow(5), border=NA,
        xlab="Individuals", ylab="Ancestry", main="(c) K = 5",
        cex.main=1.2, cex.lab=1.1)

dev.off()
cat("Figure 4 (ADMIXTURE multi-panel) saved\n")

# ============================================================
# FIGURE 5: DAPC multi-panel
# K=2 scatter (a) + K=5 scatter (b) + K=5 membership (c)
# Use manual plotting to avoid scatter.dapc layout conflicts
# ============================================================
pdf("Figure5_dapc.pdf", width=10, height=8)
layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE), heights=c(1, 0.6))

# Panel (a): DAPC K=2 — manual scatter
par(mar=c(4.5,4.5,2,2))
col2 <- c("tomato","steelblue")
pca_df2 <- as.data.frame(pca_x[, 1:2])
plot(dapc2$ind.coord[,1], rep(0, N), pch=19, cex=1.2,
     col=col2[grp2$grp], xlab="Discriminant function 1", ylab="",
     main="(a) DAPC K = 2", cex.main=1.3, yaxt="n")
legend("bottomright", legend=paste("Cluster", 1:2), fill=col2, cex=0.9)
box(lwd=1)

# Panel (b): DAPC K=5 — manual scatter
par(mar=c(4.5,4.5,2,2))
col5 <- rainbow(5)
plot(dapc5$ind.coord[,1], dapc5$ind.coord[,2], pch=19, cex=1.2,
     col=col5[grp5$grp],
     xlab=paste0("Discriminant function 1 (", round(dapc5$eig[1]/sum(dapc5$eig)*100, 1), "%)"),
     ylab=paste0("Discriminant function 2 (", round(dapc5$eig[2]/sum(dapc5$eig)*100, 1), "%)"),
     main="(b) DAPC K = 5", cex.main=1.3)
legend("bottomright", legend=paste("Cluster", 1:5), fill=col5, cex=0.8)
box(lwd=1)

# Panel (c): K=5 membership
par(mar=c(4,4,2,2))
matplot(dapc5$posterior, type="l", lty=1, lwd=1.5,
        col=col5,
        xlab="Individuals", ylab="Membership probability",
        main="(c) Membership probabilities (K = 5)", cex.main=1.3)

dev.off()
cat("Figure 5 (DAPC multi-panel) saved\n")

# ============================================================
# FIGURE 6: Kinship heatmap (single panel, already exists)
# Just verify it's using the 379-sample matrix
# ============================================================
kins_raw <- read.table("rice_kinship_379.rel", fill=TRUE, col.names=1:379)
kins <- as.matrix(kins_raw)
for (i in 1:nrow(kins)) {
  for (j in i:ncol(kins)) {
    if (i == j) next
    kins[i, j] <- kins[j, i]
  }
}

pdf("Figure6_kinship.pdf", width=10, height=9)
heatmap(kins, col=colorRampPalette(c("white","yellow","red"))(100),
        main="Kinship Matrix Heatmap (379 accessions, 26,474 SNPs)",
        cexRow=0.4, cexCol=0.4)
dev.off()
cat("Figure 6 (Kinship heatmap) saved\n")

# ============================================================
# SUPPLEMENTARY FIGURE S1: ADMIXTURE K=12
# ============================================================
q12 <- as.matrix(read.table("rice_pruned3.12.Q"))
rownames(q12) <- samples
max_pop12 <- apply(q12, 1, which.max)
ord12 <- order(max_pop12)

pdf("Supplementary_Figure_S1.pdf", width=12, height=4)
par(mar=c(4,4,1,1))
barplot(t(q12[ord12, ]), col=rainbow(12), border=NA,
        xlab="Individuals", ylab="Ancestry",
        main="ADMIXTURE K = 12 (Supplementary)")
dev.off()
cat("Supplementary Figure S1 (ADMIXTURE K=12) saved\n")

cat("\nAll multi-panel figures and supplementary figure generated.\n")
