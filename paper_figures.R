# ============================================================
# Publication-quality figures for the methodological paper
# ============================================================
library(adegenet)
library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)

workdir <- "/home/yehia/sensitivity"
setwd(workdir)

# ============================================================
# FIGURE 1: Workflow Diagram
# ============================================================
pdf("/mnt/d/project genomic analysis/Figure_1_workflow.pdf", width=10, height=11)
par(mar=c(0,0,0,0), bg="white")
plot(0, type="n", xlim=c(0,10), ylim=c(0,14), axes=FALSE, xlab="", ylab="")

# Helper: draw rounded box
draw_box <- function(x, y, w, h, col="white", text_col="black", cex=0.8, label="") {
  rect(x-w/2, y-h/2, x+w/2, y+h/2, col=col, border="gray40", lwd=1.5)
  text(x, y, label, cex=cex, col=text_col, font=2)
}

# Helper: draw arrow
draw_arrow <- function(x1, y1, x2, y2) {
  arrows(x1, y1, x2, y2, length=0.08, lwd=1.5, col="gray40")
}

# Title
text(5, 13.5, "Reproducible Population-Genomic Analysis Framework", 
     cex=1.3, font=2, col="navy")

# Step 1: Raw Dataset
draw_box(5, 12, 5, 0.6, col="#E8F5E9", cex=0.8, 
         label="Raw RDP1 Dataset\n413 accessions, 26,474 SNPs")
draw_arrow(5, 11.7, 5, 10.8)

# Step 2: QC Sensitivity
draw_box(5, 10, 7, 1.0, col="#FFF3E0", cex=0.9,
         label="QC Sensitivity Analysis")
# QC sub-boxes
draw_box(2, 9.2, 2.5, 0.5, col="#FFE0B2", cex=0.6,
         label="Mind\n(none, 0.02, 0.05, 0.10)")
draw_box(4.5, 9.2, 2.5, 0.5, col="#FFE0B2", cex=0.6,
         label="MAF\n(none, 0.01, 0.03, 0.05, 0.10)")
draw_box(7, 9.2, 2.5, 0.5, col="#FFE0B2", cex=0.6,
         label="Geno\n(none, 0.02, 0.05, 0.10, 0.20)")
draw_box(9.5, 9.2, 2.5, 0.5, col="#FFE0B2", cex=0.6,
         label="LD Pruning\n(none, 0.8, 0.5, 0.2, alt)")
text(5, 8.6, "16 scenarios varying one parameter per dimension", cex=0.6, col="gray50")

draw_arrow(5, 8.1, 5, 7.2)

# Step 3: Population Genomic Analyses
draw_box(5, 6.5, 7, 0.6, col="#E3F2FD", cex=0.9,
         label="Population Genomic Analyses (across all 16 scenarios)")

# Sub-analysis boxes
draw_box(1.5, 5.5, 2, 0.5, col="#BBDEFB", cex=0.55, label="PCA\n(10 PCs)")
draw_box(4, 5.5, 2, 0.5, col="#BBDEFB", cex=0.55, label="ADMIXTURE\nK = 1-10")
draw_box(6.5, 5.5, 2, 0.5, col="#BBDEFB", cex=0.55, label="FST / AMOVA\n5 clusters")
draw_box(9, 5.5, 2, 0.5, col="#BBDEFB", cex=0.55, label="Diversity\nPIC, He, Ho")
draw_box(2.5, 4.7, 2, 0.5, col="#BBDEFB", cex=0.55, label="DAPC\ncross-validated")
draw_box(5.5, 4.7, 2, 0.5, col="#BBDEFB", cex=0.55, label="LD Decay\nr² vs distance")
draw_box(8.5, 4.7, 2.5, 0.5, col="#BBDEFB", cex=0.55, label="Cluster Stability\naligned labels")

draw_arrow(5, 4.2, 5, 3.4)

# Step 4: Sensitivity Comparison
draw_box(5, 2.8, 7, 0.5, col="#F3E5F5", cex=0.85,
         label="Sensitivity Comparison across 16 Scenarios")
# Metrics
draw_box(1.5, 2.1, 2.5, 0.35, col="#E1BEE7", cex=0.5,
         label="Confusion Matrices")
draw_box(4.2, 2.1, 2.5, 0.35, col="#E1BEE7", cex=0.5,
         label="ARI / NMI")
draw_box(6.8, 2.1, 2.5, 0.35, col="#E1BEE7", cex=0.5,
         label="PC Correlations")
draw_box(9.3, 2.1, 1.5, 0.35, col="#E1BEE7", cex=0.5,
         label="FST / AMOVA")

draw_arrow(5, 1.7, 5, 1.1)

# Step 5: Recommendations
draw_box(5, 0.6, 6, 0.45, col="#E8F5E9", cex=0.75,
         label="Best-Practice Recommendations + Reproducible Workflow")

# Side annotation
text(0.3, 7, "Performed\nfor each\nscenario", cex=0.5, col="gray50", adj=0)

dev.off()
cat("Figure 1 (workflow) saved\n")

# ============================================================
# TABLE S15: Best-Practice Recommendations
# ============================================================
cat("\nGenerating recommendation table...\n")

# Load sensitivity master table
master <- tryCatch(read.csv("/mnt/d/project genomic analysis/Table_S4_sensitivity_master.csv",
                            stringsAsFactors=FALSE), error=function(e) NULL)
# Load diversity
diversity <- read.csv("/mnt/d/project genomic analysis/Table_S10_genetic_diversity.csv",
                      stringsAsFactors=FALSE)
# Load FST
fst <- read.csv("/mnt/d/project genomic analysis/Table_S11_FST_comparison.csv",
                stringsAsFactors=FALSE)
# Load ARI
ari <- read.csv("/mnt/d/project genomic analysis/Table_S7_pairwise_ARI.csv",
                stringsAsFactors=FALSE, row.names=1)
# Load LD decay
ld <- read.csv("/mnt/d/project genomic analysis/Table_S13_LD_decay.csv",
               stringsAsFactors=FALSE)

# Compute deviations from baseline
baseline_div <- diversity[diversity$Short == "S3", ]

rec_table <- data.frame(
  Dimension = character(),
  Threshold = character(),
  Impact_on_SNP_retention = character(),
  Impact_on_Sample_retention = character(),
  Impact_on_Diversity = character(),
  Impact_on_FST = character(),
  Impact_on_Cluster_ARI = character(),
  Impact_on_LD = character(),
  Recommendation = character(),
  stringsAsFactors=FALSE
)

# Sample filtering
rec_table <- rbind(rec_table, data.frame(
  Dimension = "Sample filtering", Threshold = "mind = none",
  SNP_retention = "1,166 SNPs", Sample_retention = "413 (100%)",
  Diversity = "Baseline", FST = "0.414",
  Cluster_ARI = "0.995", LD = "Not applicable (no SNP filter)",
  Recommendation = "Not recommended – retains noisy samples",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "Sample filtering", Threshold = "mind = 0.10 (after QC)",
  SNP_retention = "1,161 SNPs", Sample_retention = "409 (99%)",
  Diversity = "Near-identical", FST = "0.416",
  Cluster_ARI = "0.995", LD = "Not applicable",
  Recommendation = "Acceptable – minimal sample loss",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "Sample filtering", Threshold = "mind = 0.10 (baseline)",
  SNP_retention = "1,187 SNPs", Sample_retention = "379 (92%)",
  Diversity = "Reference", FST = "0.423",
  Cluster_ARI = "1.000", LD = "Reference",
  Recommendation = "Recommended – balances retention and quality",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "Sample filtering", Threshold = "mind = 0.05",
  SNP_retention = "1,064 SNPs", Sample_retention = "293 (71%)",
  Diversity = "Similar", FST = "0.383 (-9.4%)",
  Cluster_ARI = "0.580", LD = "Not applicable",
  Recommendation = "Too stringent – loses 29% of samples, alters FST",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "Sample filtering", Threshold = "mind = 0.02",
  SNP_retention = "597 SNPs", Sample_retention = "98 (24%)",
  Diversity = "Reduced (fewer SNPs)", FST = "0.293 (-30.7%)",
  Cluster_ARI = "0.561", LD = "Not applicable",
  Recommendation = "Severely biased – do not use",
  stringsAsFactors=FALSE
))

# MAF filtering
rec_table <- rbind(rec_table, data.frame(
  Dimension = "MAF filtering", Threshold = "MAF = none",
  SNP_retention = "1,638 SNPs (+38%)", Sample_retention = "383 (100%)",
  Diversity = "Includes rare variants", FST = "0.348 (-17.8%)",
  Cluster_ARI = "0.995", LD = "Not applicable",
  Recommendation = "Rare variants inflate noise – not recommended",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "MAF filtering", Threshold = "MAF = 0.01",
  SNP_retention = "1,601 SNPs (+35%)", Sample_retention = "383 (100%)",
  Diversity = "Slightly inflated", FST = "0.353 (-16.7%)",
  Cluster_ARI = "0.992", LD = "Not applicable",
  Recommendation = "Marginal – rare variants dilute structure",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "MAF filtering", Threshold = "MAF = 0.03",
  SNP_retention = "1,386 SNPs (+17%)", Sample_retention = "383 (100%)",
  Diversity = "Near-baseline", FST = "0.394 (-6.9%)",
  Cluster_ARI = "0.997", LD = "Not applicable",
  Recommendation = "Acceptable – modest SNP gain, FST stable",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "MAF filtering", Threshold = "MAF = 0.05 (baseline)",
  SNP_retention = "1,187 SNPs", Sample_retention = "379 (100%)",
  Diversity = "Reference", FST = "0.423",
  Cluster_ARI = "1.000", LD = "Reference",
  Recommendation = "Recommended – standard for population structure",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "MAF filtering", Threshold = "MAF = 0.10",
  SNP_retention = "844 SNPs (-29%)", Sample_retention = "383 (100%)",
  Diversity = "Reduced PIC, Ho inflated", FST = "0.449 (+6.1%)",
  Cluster_ARI = "0.863", LD = "Not applicable",
  Recommendation = "Too stringent for structure – loses informative SNPs",
  stringsAsFactors=FALSE
))

# GENO filtering
rec_table <- rbind(rec_table, data.frame(
  Dimension = "GENO filtering", Threshold = "GENO = none",
  SNP_retention = "1,437 SNPs (+21%)", Sample_retention = "383 (100%)",
  Diversity = "Slightly inflated", FST = "0.415 (-1.9%)",
  Cluster_ARI = "0.995", LD = "Not applicable",
  Recommendation = "Acceptable – small effect on structure",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "GENO filtering", Threshold = "GENO = 0.20",
  SNP_retention = "1,404 SNPs (+18%)", Sample_retention = "383 (100%)",
  Diversity = "Similar", FST = "Not computed",
  Cluster_ARI = "Not computed", LD = "Not applicable",
  Recommendation = "Acceptable – retains more SNPs",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "GENO filtering", Threshold = "GENO = 0.10",
  SNP_retention = "1,299 SNPs (+9%)", Sample_retention = "383 (100%)",
  Diversity = "Similar", FST = "0.397 (-6.2%)",
  Cluster_ARI = "0.858", LD = "Not applicable",
  Recommendation = "Acceptable – modest SNP gain",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "GENO filtering", Threshold = "GENO = 0.05 (baseline)",
  SNP_retention = "1,187 SNPs", Sample_retention = "379 (100%)",
  Diversity = "Reference", FST = "0.423",
  Cluster_ARI = "1.000", LD = "Reference",
  Recommendation = "Recommended standard",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "GENO filtering", Threshold = "GENO = 0.02",
  SNP_retention = "982 SNPs (-17%)", Sample_retention = "383 (100%)",
  Diversity = "Reduced", FST = "0.419 (-1.0%)",
  Cluster_ARI = "0.992", LD = "Not applicable",
  Recommendation = "Stringent but FST stable",
  stringsAsFactors=FALSE
))

# LD pruning
rec_table <- rbind(rec_table, data.frame(
  Dimension = "LD pruning", Threshold = "none (all SNPs)",
  SNP_retention = "26,474 SNPs", Sample_retention = "383 (100%)",
  Diversity = "Near-identical", FST = "0.609 (+44.0%)",
  Cluster_ARI = "0.989", LD = "r² = 0.43 at 0–10kb (48% > 0.2)",
  Recommendation = "Not recommended – inflates FST, violates independence",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "LD pruning", Threshold = "r² < 0.8",
  SNP_retention = "9,609 SNPs", Sample_retention = "383 (100%)",
  Diversity = "Near-identical", FST = "0.497 (+17.5%)",
  Cluster_ARI = "0.989", LD = "r² = 0.21 at 0–10kb (23% > 0.2)",
  Recommendation = "Insufficient pruning – FST still inflated",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "LD pruning", Threshold = "r² < 0.5",
  SNP_retention = "4,278 SNPs", Sample_retention = "383 (100%)",
  Diversity = "Near-identical", FST = "0.434 (+2.6%)",
  Cluster_ARI = "0.995", LD = "r² = 0.14 at 0–10kb (13% > 0.2)",
  Recommendation = "Acceptable – FST near-baseline",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "LD pruning", Threshold = "50/5/0.2 (baseline)",
  SNP_retention = "1,187 SNPs", Sample_retention = "379 (100%)",
  Diversity = "Reference", FST = "0.423 (reference)",
  Cluster_ARI = "1.000", LD = "r² = 0.07 at 0–10kb (5.8% > 0.2)",
  Recommendation = "Recommended – standard for population structure",
  stringsAsFactors=FALSE
))
rec_table <- rbind(rec_table, data.frame(
  Dimension = "LD pruning", Threshold = "100/10/0.2",
  SNP_retention = "822 SNPs (-31%)", Sample_retention = "383 (100%)",
  Diversity = "Similar", FST = "0.353 (-16.6%)",
  Cluster_ARI = "0.799", LD = "r² = 0.07 at 0–10kb (1.4% > 0.2)",
  Recommendation = "Overly aggressive – loses SNPs, reduces FST",
  stringsAsFactors=FALSE
))

write.csv(rec_table, "/mnt/d/project genomic analysis/Table_S15_recommendations.csv", row.names=FALSE)
cat("Table S15 (recommendations) written:", nrow(rec_table), "rows\n")

# ============================================================
# FIGURE 16: Summary of Recommendations (visual dashboard)
# ============================================================
pdf("/mnt/d/project genomic analysis/Figure_16_recommendations.pdf", width=10, height=7)
par(mfrow=c(2,2), mar=c(5,4,3,1))

# 1. SNP retention by dimension
dims <- c("Sample", "Sample", "Sample", "Sample", "Sample",
          "MAF", "MAF", "MAF", "MAF", "MAF",
          "GENO", "GENO", "GENO", "GENO", "GENO",
          "LD", "LD", "LD", "LD", "LD")
thresholds <- rec_table$Threshold
snps <- c(1166, 1161, 1187, 1064, 597,
          1638, 1601, 1386, 1187, 844,
          1437, 1404, 1299, 1187, 982,
          26474, 9609, 4278, 1187, 822)
barplot(snps, names.arg=thresholds, las=2, cex.names=0.5, col=rep(c("skyblue","salmon","lightgreen","plum"), each=5),
        main="SNP Retention by QC Threshold", ylab="Number of SNPs")

# 2. FST change  
fst_vals <- c(0.414, 0.416, 0.423, 0.383, 0.293,
              0.348, 0.353, 0.394, 0.423, 0.449,
              NA, NA, 0.397, 0.423, 0.419,
              0.609, 0.497, 0.434, 0.423, 0.353)
barplot(fst_vals, names.arg=thresholds, las=2, cex.names=0.5,
        col=rep(c("skyblue","salmon","lightgreen","plum"), each=5),
        main="Mean FST by QC Threshold", ylab="FST")
abline(h=0.423, lty=2, col="red", lwd=2)

# 3. Cluster stability (ARI)
ari_vals <- c(0.995, 0.995, 1.0, 0.580, 0.561,
              0.995, 0.992, 0.997, 1.0, 0.863,
              0.995, NA, 0.858, 1.0, 0.992,
              0.989, 0.989, 0.995, 1.0, 0.799)
barplot(ari_vals, names.arg=thresholds, las=2, cex.names=0.5,
        col=rep(c("skyblue","salmon","lightgreen","plum"), each=5),
        main="Cluster Stability (ARI vs Baseline)", ylab="ARI", ylim=c(0,1))
abline(h=0.9, lty=2, col="red", lwd=2)

# 4. LD at short range
ld_r2 <- c(NA, NA, NA, NA, NA,
           NA, NA, NA, NA, NA,
           NA, NA, NA, NA, NA,
           0.43, 0.21, 0.14, 0.074, 0.070)
barplot(ld_r2, names.arg=thresholds, las=2, cex.names=0.5,
        col=rep("plum", 20),
        main="LD Decay: Mean r² at 0-10kb", ylab="Mean r²")
abline(h=0.074, lty=2, col="red", lwd=2)

dev.off()
cat("Figure 16 (recommendations dashboard) saved\n")

cat("\n===== ALL PAPER FIGURES GENERATED =====\n")
