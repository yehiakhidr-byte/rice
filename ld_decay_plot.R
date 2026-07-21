library(ggplot2)

workdir <- "/home/yehia/sensitivity"
setwd(workdir)

cat("Processing LD decay data...\n")

ld_scenarios <- c("sens_D4_LDnone", "sens_D4_LD08", "sens_D4_LD05", "rice_pruned3", "sens_D4_LD100")
ld_labels <- c("L0:LDnone", "L1:LD08", "L2:LD05", "L3:LD02(baseline)", "L4:LD100_10_02")
ld_n_snps <- c(26474, 9609, 4278, 1187, 822)

# Distance bins (kb)
bins <- c(0, 10, 50, 100, 200, 500, 1000, 5000)
bin_labels <- c("0-10", "10-50", "50-100", "100-200", "200-500", "500-1000", "1000-5000")

ld_decay_results <- data.frame()

for (i in seq_along(ld_scenarios)) {
  f <- ld_scenarios[i]
  lbl <- ld_labels[i]
  n_snp <- ld_n_snps[i]
  
  ld_file <- paste0(f, "_ld_decay.ld")
  if (!file.exists(ld_file)) {
    cat("  Missing:", ld_file, "\n")
    next
  }
  
  cat("  Reading", ld_file, "...\n")
  ld <- read.table(ld_file, header=TRUE, stringsAsFactors=FALSE)
  
  # Compute distance in kb
  ld$dist_kb <- abs(ld$BP_B - ld$BP_A) / 1000
  
  # Bin by distance
  ld$bin <- cut(ld$dist_kb, breaks=bins, labels=bin_labels, include.lowest=TRUE)
  
  # Mean r² per bin
  mean_r2 <- tapply(ld$R2, ld$bin, mean, na.rm=TRUE)
  n_pairs <- tapply(ld$R2, ld$bin, length)
  
  result <- data.frame(
    Scenario = lbl,
    n_SNPs = n_snp,
    Bin = names(mean_r2),
    Midpoint_kb = c(5, 30, 75, 150, 350, 750, 3000),
    Mean_r2 = as.numeric(mean_r2),
    N_pairs = as.numeric(n_pairs),
    stringsAsFactors=FALSE
  )
  ld_decay_results <- rbind(ld_decay_results, result)
  
  cat("    Total pairs:", nrow(ld), "\n")
}

# Write table
write.csv(ld_decay_results, "Table_S13_LD_decay.csv", row.names=FALSE)

# Figure 14: LD decay plot
cat("Generating Figure 14 (LD decay)...\n")

pdf("Figure_14_LD_decay.pdf", width=8, height=6)

# Color by scenario
scenario_colors <- c("L0:LDnone"="#E41A1C", "L1:LD08"="#377EB8", 
                     "L2:LD05"="#4DAF4A", "L3:LD02(baseline)"="#984EA3",
                     "L4:LD100_10_02"="#FF7F00")

plot(ld_decay_results$Midpoint_kb, ld_decay_results$Mean_r2, type="n",
     xlab="Distance (kb)", ylab=expression(Mean~r^2),
     main="LD Decay Across Pruning Scenarios",
     log="x", xlim=c(3, 5000), ylim=c(0, max(ld_decay_results$Mean_r2, na.rm=TRUE) * 1.1),
     xaxt="n")

# Custom x-axis
axis(1, at=c(5, 10, 50, 100, 200, 500, 1000, 3000), 
     labels=c("5", "10", "50", "100", "200", "500", "1000", "3000"))

for (sc in unique(ld_decay_results$Scenario)) {
  sub <- ld_decay_results[ld_decay_results$Scenario == sc, ]
  sub <- sub[order(sub$Midpoint_kb), ]
  lines(sub$Midpoint_kb, sub$Mean_r2, type="o", pch=19, 
        col=scenario_colors[sc], lwd=2, cex=1.2)
}

legend("topright", legend=names(scenario_colors), col=scenario_colors,
       lwd=2, pch=19, cex=0.8, title="Scenario (n SNPs)")

# Add inset table with n_SNPs
legend("bottomleft", 
       legend=paste0(ld_labels, " (", ld_n_snps, " SNPs)"),
       cex=0.6, bty="n")

dev.off()
cat("Figure 14 (LD decay) saved\n")

# Also add LD decay line plot to the standard output PDFs
# Comparison of LD at short distances
cat("Generating LD decay summary...\n")

# Key metric: mean r² at short distances (0-10kb)
short_range <- ld_decay_results[ld_decay_results$Bin == "0-10", ]
cat("Mean r² at 0-10kb:\n")
print(short_range[, c("Scenario", "Mean_r2", "n_SNPs")])

# Additional comparison: how much LD is removed by pruning
cat("\nProportion of pairs with r² > 0.2:\n")
for (i in seq_along(ld_scenarios)) {
  f <- ld_scenarios[i]
  lbl <- ld_labels[i]
  ld_file <- paste0(f, "_ld_decay.ld")
  ld <- read.table(ld_file, header=TRUE, stringsAsFactors=FALSE)
  pct_high_ld <- mean(ld$R2 > 0.2, na.rm=TRUE) * 100
  cat(sprintf("  %s: %.1f%%\n", lbl, pct_high_ld))
}

cat("\n===== LD DECAY COMPLETE =====\n")
