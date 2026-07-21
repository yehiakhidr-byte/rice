# Kinship summary for 379-sample matrix
kins_raw <- read.table("rice_kinship_379.rel", fill=TRUE, col.names=1:379)
kins <- as.matrix(kins_raw)
cat("Kinship matrix: 379 x 379\n")

# Matrix should be lower-triangular with diagonal; symmetrize
for (i in 1:nrow(kins)) {
  for (j in i:ncol(kins)) {
    if (i == j) next
    kins[i, j] <- kins[j, i]
  }
}
upper <- kins[upper.tri(kins)]

cat("=== Table: Kinship Summary (379 samples) ===\n")
cat(sprintf("%-30s %s\n", "Statistic", "Value"))
cat(strrep("-", 45), "\n")
cat(sprintf("%-30s %d x %d\n", "Matrix dimensions", nrow(kins), ncol(kins)))
cat(sprintf("%-30s %.4f\n", "Mean", mean(upper)))
cat(sprintf("%-30s %.4f\n", "SD", sd(upper)))
cat(sprintf("%-30s %.4f\n", "Min", min(upper)))
cat(sprintf("%-30s %.4f\n", "Max", max(upper)))
cat(sprintf("%-30s %.4f\n", "Median", median(upper)))
cat(sprintf("%-30s %.2f%%\n", "Negative values", sum(upper < 0) / length(upper) * 100))
cat(sprintf("%-30s %.2f%%\n", ">0.25 (close kin)", sum(upper > 0.25) / length(upper) * 100))
cat(sprintf("%-30s %.2f%%\n", ">0.125 (2nd degree)", sum(upper > 0.125) / length(upper) * 100))

# Regenerate plots
library(RColorBrewer)
pdf("kinship_heatmap_v4.pdf", width=10, height=9)
heatmap(kins, col=colorRampPalette(c("white","yellow","red"))(100),
        main="Kinship Matrix Heatmap (379 samples, 26,474 SNPs)", cexRow=0.5, cexCol=0.5)
dev.off()
cat("Kinship heatmap (379) saved\n")

pdf("kinship_distribution_v4.pdf", width=6, height=4)
hist(upper, breaks=50, col="steelblue",
     main="Pairwise Kinship Coefficients (379 samples)", xlab="Kinship coefficient")
dev.off()
cat("Kinship distribution (379) saved\n")
