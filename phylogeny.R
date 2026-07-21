library(adegenet)
library(ape)

setwd("/mnt/d/project genomic analysis")

# Load genotype data
obj <- read.PLINK("rice_pruned3.raw", parallel=FALSE)
gt <- tab(obj)
for(i in 1:ncol(gt)){gt[is.na(gt[,i]), i] <- mean(gt[,i], na.rm=TRUE)}
fam <- read.table("rice_pruned3.fam")
samples <- fam[, 2]

# Compute pairwise genetic distance
cat("Computing genetic distances...\n")
dist_mat <- dist(gt, method = "euclidean")

# Build NJ tree
cat("Building Neighbor-Joining tree...\n")
nj_tree <- nj(dist_mat)

# Load ADMIXTURE K=5 assignments for coloring
q5 <- as.matrix(read.table("rice_pruned3.5.Q"))
admix_clust <- apply(q5, 1, which.max)

clust_colors <- c("red", "blue", "green", "gold", "purple")
tip_colors <- clust_colors[admix_clust]

cat("Cluster sizes (ADMIXTURE K=5):\n")
print(table(admix_clust))

# Unrooted tree (main figure)
pdf("Figure7_phylogeny.pdf", width=12, height=14)
par(mar=c(1,1,1,1))
plot(nj_tree, type="unrooted", tip.color=tip_colors, cex=0.5,
     main="Neighbor-Joining Phylogenetic Tree of 379 Rice Accessions")
legend("bottomleft", legend=paste("Cluster", 1:5), fill=clust_colors,
       cex=0.8, title="ADMIXTURE K=5")
add.scale.bar(cex=0.7)
dev.off()
cat("Unrooted tree saved: Figure7_phylogeny.pdf\n")

# Ladderized phylogram (supplementary)
nj_lad <- ladderize(nj_tree)

pdf("Supplementary_Figure_S2_phylogeny.pdf", width=10, height=12)
par(mar=c(1,1,1,1))
plot(nj_lad, type="phylogram", tip.color=tip_colors, cex=0.4,
     direction="rightwards",
     main="Neighbor-Joining Tree of 379 Rice Accessions")
legend("topleft", legend=paste("Cluster", 1:5), fill=clust_colors,
       cex=0.8, title="ADMIXTURE K=5")
dev.off()
cat("Ladderized tree saved: Supplementary_Figure_S2_phylogeny.pdf\n")
