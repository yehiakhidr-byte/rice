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
n_pcs <- min(40, N - 1)
pca_x <- pca$x[, 1:n_pcs]

# ========== ADMIXTURE K=2 BAR PLOT ==========
q2 <- as.matrix(read.table("rice_pruned3.2.Q"))
rownames(q2) <- samples
max_pop2 <- apply(q2, 1, which.max)
ord2 <- order(max_pop2)
q2_ord <- q2[ord2, ]

pdf("admixture_K2_barplot_v4.pdf", width=12, height=4)
par(mar=c(4,4,1,1))
barplot(t(q2_ord), col=c("tomato","steelblue"), border=NA,
        xlab="Individuals", ylab="Ancestry",
        main="ADMIXTURE Ancestry Proportions (K=2, 379 samples)")
dev.off()
cat("ADMIXTURE K=2 bar plot saved\n")

# ========== DAPC K=2 (forced) ==========
set.seed(42)
grp2 <- find.clusters(pca_x, max.n.clust=20, n.iter=1e6, n.start=10,
                      stat="BIC", n.pca=n_pcs, n.clust=2)
cat("DAPC K=2 cluster sizes:\n")
print(table(grp2$grp))

dapc2 <- dapc(pca_x, grp2$grp, n.pca=min(30, N-1, ncol(pca_x)), n.da=1)

pdf("dapc_scatter_K2_v4.pdf", width=7, height=6)
scatter(dapc2, scree.da=FALSE, legend=TRUE,
        posi.leg="bottomright", clab=0, cex=1.5,
        main="DAPC K=2 (379 samples)")
dev.off()
cat("DAPC K=2 scatter saved\n")

pdf("dapc_membership_K2_v4.pdf", width=10, height=4)
matplot(dapc2$posterior, type="l", lty=1, lwd=1,
        col=c("tomato","steelblue"),
        xlab="Individuals", ylab="Membership",
        main="DAPC Membership Probabilities (K=2)")
dev.off()
cat("DAPC K=2 membership saved\n")

# ========== DAPC K=5 (already in main script but rerun for consistency) ==========
set.seed(42)
grp5 <- find.clusters(pca_x, max.n.clust=20, n.iter=1e6, n.start=10,
                      stat="BIC", n.pca=n_pcs, n.clust=5)
cat("\nDAPC K=5 cluster sizes:\n")
print(table(grp5$grp))

dapc5 <- dapc(pca_x, grp5$grp, n.pca=min(30, N-1, ncol(pca_x)), n.da=4)

# ========== CORRESPONDENCE: ADMIXTURE K=2 vs DAPC K=2 ==========
# Cross-tabulate admixture K=2 assignment (max Q) vs DAPC K=2 cluster
admix_k2_clust <- apply(q2, 1, which.max)
cat("\nCross-tabulation: ADMIXTURE K=2 vs DAPC K=2\n")
print(table(ADMIXTURE=admix_k2_clust, DAPC=grp2$grp))

# ========== CORRESPONDENCE: ADMIXTURE K=5 vs DAPC K=5 ==========
q5 <- as.matrix(read.table("rice_pruned3.5.Q"))
admix_k5_clust <- apply(q5, 1, which.max)
cat("\nCross-tabulation: ADMIXTURE K=5 vs DAPC K=5\n")
print(table(ADMIXTURE=admix_k5_clust, DAPC=grp5$grp))

cat("\nAll K=2 and correspondence analyses complete.\n")
