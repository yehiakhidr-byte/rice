library(adegenet)

setwd("/mnt/d/project genomic analysis")

obj <- read.PLINK("rice_pruned3.raw", parallel=FALSE)
gt <- tab(obj)

cat("=== MARKER DIVERSITY SUMMARY ===\n")
cat(sprintf("Total SNPs: %d\n", ncol(gt)))
cat(sprintf("Total individuals: %d\n\n", nrow(gt)))

# Minor allele frequency
maf <- apply(gt, 2, function(x) {
  f <- mean(x, na.rm=TRUE) / 2
  min(f, 1 - f)
})
cat(sprintf("Mean MAF: %.4f\n", mean(maf, na.rm=TRUE)))
cat(sprintf("Median MAF: %.4f\n", median(maf, na.rm=TRUE)))
cat(sprintf("MAF < 0.05: %d (%.1f%%)\n", sum(maf < 0.05, na.rm=TRUE),
            sum(maf < 0.05, na.rm=TRUE) / length(maf) * 100))
cat(sprintf("MAF 0.05–0.10: %d (%.1f%%)\n", sum(maf >= 0.05 & maf < 0.10, na.rm=TRUE),
            sum(maf >= 0.05 & maf < 0.10, na.rm=TRUE) / length(maf) * 100))
cat(sprintf("MAF 0.10–0.30: %d (%.1f%%)\n", sum(maf >= 0.10 & maf < 0.30, na.rm=TRUE),
            sum(maf >= 0.10 & maf < 0.30, na.rm=TRUE) / length(maf) * 100))
cat(sprintf("MAF >= 0.30: %d (%.1f%%)\n\n", sum(maf >= 0.30, na.rm=TRUE),
            sum(maf >= 0.30, na.rm=TRUE) / length(maf) * 100))

# Observed and expected heterozygosity
ho <- apply(gt, 2, function(x) mean(x == 1, na.rm=TRUE))
he <- apply(gt, 2, function(x) {
  p <- mean(x, na.rm=TRUE) / 2
  2 * p * (1 - p)
})
cat(sprintf("Mean observed heterozygosity (Ho): %.4f\n", mean(ho, na.rm=TRUE)))
cat(sprintf("SD Ho: %.4f\n", sd(ho, na.rm=TRUE)))
cat(sprintf("Mean expected heterozygosity (He): %.4f\n", mean(he, na.rm=TRUE)))
cat(sprintf("SD He: %.4f\n\n", sd(he, na.rm=TRUE)))

# PIC (Polymorphic Information Content) = He for biallelic markers
cat(sprintf("Mean PIC: %.4f\n\n", mean(he, na.rm=TRUE)))

# Missing rate
miss <- apply(gt, 2, function(x) sum(is.na(x)) / length(x))
cat(sprintf("Mean missing rate per SNP: %.4f\n", mean(miss)))
cat(sprintf("SNPs with 0%% missing: %d (%.1f%%)\n", sum(miss == 0),
            sum(miss == 0) / length(miss) * 100))
cat(sprintf("SNPs with >5%% missing: %d (%.1f%%)\n\n", sum(miss > 0.05),
            sum(miss > 0.05) / length(miss) * 100))

# By chromosome
cat("=== DIVERSITY BY CHROMOSOME ===\n")
bim <- read.table("rice_pruned3.bim", stringsAsFactors=FALSE)
names(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
chr_stats <- aggregate(he ~ bim$CHR, FUN=function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
chr_order <- order(as.numeric(as.character(chr_stats[,1])))
cat(sprintf("%-4s %8s %8s %8s\n", "Chr", "N_SNPs", "Mean_He", "SD_He"))
for (i in chr_order) {
  cat(sprintf("%-4s %8.0f %8.4f %8.4f\n", chr_stats[i,1], chr_stats[i,2][3],
              chr_stats[i,2][1], chr_stats[i,2][2]))
}

# MAF distribution plot
pdf("marker_diversity_maf.pdf", width=8, height=5)
par(mar=c(4.5,4.5,2,1))
hist(maf, breaks=50, col="steelblue", border="white",
     xlab="Minor allele frequency", ylab="Number of SNPs",
     main="MAF Distribution (26,474 QC-passing SNPs)")
dev.off()
cat("\nMAF distribution plot saved: marker_diversity_maf.pdf\n")

# He distribution plot
pdf("marker_diversity_he.pdf", width=8, height=5)
par(mar=c(4.5,4.5,2,1))
hist(he, breaks=50, col="coral", border="white",
     xlab="Expected heterozygosity (He)", ylab="Number of SNPs",
     main="He Distribution (26,474 QC-passing SNPs)")
dev.off()
cat("He distribution plot saved: marker_diversity_he.pdf\n")

cat("\nAll marker diversity analyses complete.\n")
