library(adegenet)
library(RColorBrewer)

setwd("/mnt/d/project genomic analysis")

# ===== 1. LOAD PASSPORT DATA =====
passport <- read.csv("/mnt/d/Rice accessions/RDP1_full_crossref.csv", stringsAsFactors=FALSE)
names(passport)[1] <- "NSFTV_ID"
# Create numeric key from NSFTV_ID (e.g., "NSFTV_1" -> 1)
passport$NSFTV_num <- as.numeric(gsub("NSFTV_", "", passport$NSFTV_ID))

# ===== 2. LOAD FAM FILE (379 samples) =====
fam <- read.table("rice_pruned3.fam", stringsAsFactors=FALSE)
names(fam) <- c("Name", "NSFTV_num", "V3", "V4", "V5", "V6")

# Merge fam with passport
annot <- merge(fam, passport[, c("NSFTV_num", "Subpopulation", "Country", "Accession_Name")],
               by="NSFTV_num", all.x=TRUE, sort=FALSE)

cat("=== PASSPORT ANNOTATION SUMMARY ===\n")
cat(sprintf("Total samples: %d\n", nrow(annot)))
cat(sprintf("Annotated: %d\n", sum(!is.na(annot$Subpopulation))))
cat(sprintf("Missing annotation: %d\n\n", sum(is.na(annot$Subpopulation))))

cat("Subpopulation counts:\n")
print(table(annot$Subpopulation))
cat("\nCountry counts (top 15):\n")
print(head(sort(table(annot$Country), decreasing=TRUE), 15))

# ===== 3. LOAD ADMIXTURE K=5 AND ALIGN LABELS =====
q5 <- as.matrix(read.table("rice_pruned3.5.Q"))
admix_clust <- apply(q5, 1, which.max)
annot$Admix_clust <- admix_clust

cat("\n=== ADMIXTURE K=5 vs PASSPORT SUBPOPULATION ===\n")
cross <- table(ADMIXTURE=annot$Admix_clust, Subpopulation=annot$Subpopulation)
print(cross)

# Find the dominant subpopulation per ADMIXTURE cluster
cat("\nDominant subpopulation per ADMIXTURE cluster:\n")
for (cl in 1:5) {
  sub <- cross[cl,]
  dominant <- names(which.max(sub))
  cat(sprintf("  Cluster %d -> %s (%d of %d)\n", cl, dominant, max(sub), sum(sub)))
}

# ===== 4. PCA COLORED BY SUBPOPULATION =====
obj <- read.PLINK("rice_pruned3.raw", parallel=FALSE)
gt <- tab(obj)
for(i in 1:ncol(gt)){gt[is.na(gt[,i]), i] <- mean(gt[,i], na.rm=TRUE)}
pca <- prcomp(gt, scale=FALSE, center=TRUE)
pve <- summary(pca)$importance[2, 1:4] * 100

subpop_colors <- c("IND"="red", "AUS"="blue", "TEJ"="green", "TRJ"="gold", "AROMATIC"="purple", "ADMIX"="gray50")
annot$Color <- subpop_colors[annot$Subpopulation]
annot$Color[is.na(annot$Color)] <- "gray"

# PCA by subpopulation
pdf("pca_by_subpop.pdf", width=8, height=7)
par(mar=c(4.5,4.5,2,1))
plot(pca$x[,1], pca$x[,2], pch=19, cex=0.7, col=annot$Color,
     xlab=paste0("PC1 (", round(pve[1], 1), "%)"),
     ylab=paste0("PC2 (", round(pve[2], 1), "%)"),
     main="PCA of 379 Rice Accessions Colored by Subpopulation")
legend("topright", legend=names(subpop_colors), fill=subpop_colors, cex=0.7)
dev.off()
cat("\nPCA by subpopulation saved: pca_by_subpop.pdf\n")

# PCA by country (top 10 countries)
top_countries <- names(sort(table(annot$Country), decreasing=TRUE)[1:10])
annot$Country_group <- ifelse(annot$Country %in% top_countries, annot$Country, "Other")
country_colors <- c(rainbow(10), "gray50")
names(country_colors) <- c(top_countries, "Other")
annot$Country_color <- country_colors[annot$Country_group]

pdf("pca_by_country.pdf", width=9, height=7)
par(mar=c(4.5,4.5,2,1))
plot(pca$x[,1], pca$x[,2], pch=19, cex=0.7, col=annot$Country_color,
     xlab=paste0("PC1 (", round(pve[1], 1), "%)"),
     ylab=paste0("PC2 (", round(pve[2], 1), "%)"),
     main="PCA of 379 Rice Accessions Colored by Country")
legend("topright", legend=names(country_colors), fill=country_colors, cex=0.6)
dev.off()
cat("PCA by country saved: pca_by_country.pdf\n")

# ===== 5. DAPC COLORED BY SUBPOPULATION =====
N <- nInd(obj)
n_pcs <- min(40, N - 1)
pca_x <- pca$x[, 1:n_pcs]

set.seed(42)
grp5 <- find.clusters(pca_x, max.n.clust=20, n.iter=1e6, n.start=10,
                      stat="BIC", n.pca=n_pcs, n.clust=5)
dapc5 <- dapc(pca_x, grp5$grp, n.pca=min(30, N-1, ncol(pca_x)), n.da=4)
annot$Dapc_clust <- grp5$grp

cat("\n=== DAPC K=5 vs PASSPORT SUBPOPULATION ===\n")
dapc_cross <- table(DAPC=annot$Dapc_clust, Subpopulation=annot$Subpopulation)
print(dapc_cross)

cat("\nDominant subpopulation per DAPC cluster:\n")
for (cl in 1:5) {
  sub <- dapc_cross[cl,]
  dominant <- names(which.max(sub))
  cat(sprintf("  Cluster %d -> %s (%d of %d)\n", cl, dominant, max(sub), sum(sub)))
}

# DAPC scatter colored by passport subpopulation
pdf("dapc_by_subpop.pdf", width=8, height=7)
par(mar=c(4.5,4.5,2,1))
plot(dapc5$ind.coord[,1], dapc5$ind.coord[,2], pch=19, cex=0.8, col=annot$Color,
     xlab=paste0("Discriminant function 1 (", round(dapc5$eig[1]/sum(dapc5$eig)*100, 1), "%)"),
     ylab=paste0("Discriminant function 2 (", round(dapc5$eig[2]/sum(dapc5$eig)*100, 1), "%)"),
     main="DAPC of 379 Rice Accessions Colored by Subpopulation")
legend("bottomright", legend=names(subpop_colors), fill=subpop_colors, cex=0.7)
dev.off()
cat("DAPC by subpopulation saved: dapc_by_subpop.pdf\n")

# ===== 6. ADMIXTURE BARPLOT WITH SUBPOPULATION LABELS =====
# Reorder by passport subpopulation
ord_subpop <- order(annot$Subpopulation, annot$Admix_clust)
q5_ord_sub <- q5[ord_subpop, ]
annot_ord <- annot[ord_subpop, ]

# Create color bar for subpopulation
subpop_bar <- subpop_colors[annot_ord$Subpopulation]
subpop_bar[is.na(subpop_bar)] <- "gray"

pdf("admixture_by_subpop.pdf", width=14, height=5)
layout(matrix(c(1,2), nrow=2), heights=c(1, 0.15))
par(mar=c(2,4,2,1))
barplot(t(q5_ord_sub), col=rainbow(5), border=NA,
        xlab="", ylab="Ancestry",
        main="ADMIXTURE K=5 Ordered by Passport Subpopulation")
par(mar=c(0,4,0,1))
image(t(1:length(subpop_bar)), col=subpop_bar, axes=FALSE)
axis(1, at=seq(0, 1, length.out=6), 
     labels=c("AROMATIC", "AUS", "ADMIX", "IND", "TEJ", "TRJ"), las=2, cex.axis=0.7)
dev.off()
cat("ADMIXTURE by subpopulation saved: admixture_by_subpop.pdf\n")

# ===== 7. SUMMARY BY COUNTRY AND SUBPOPULATION =====
cat("\n\n=== SUMMARY BY COUNTRY ===\n")
country_summary <- table(annot$Country, annot$Subpopulation)
# Show top countries
top_ctry <- names(sort(rowSums(country_summary), decreasing=TRUE)[1:10])
print(country_summary[top_ctry, ])

cat("\n\n=== SAVE ANNOTATED DATASET ===\n")
write.csv(annot[, c("NSFTV_num", "Name", "Accession_Name", "Country", "Subpopulation", "Admix_clust", "Dapc_clust")],
          "rice379_annotated.csv", row.names=FALSE)
cat("Annotated dataset saved: rice379_annotated.csv\n")
cat("\nAll passport annotation complete!\n")
