library(RColorBrewer)
library(adegenet)

setwd("/mnt/d/project genomic analysis")
SENS <- "/home/yehia/sensitivity"
BD <- "/mnt/d/Rice accessions/RiceDiversity.44K.MSU6.Genotypes_PLINK/RiceDiversity_44K_Genotypes_PLINK"

# ---- DATA ----
cv_raw <- read.table("sensitivity_cv_summary.txt", header=TRUE, stringsAsFactors=FALSE)
colnames(cv_raw) <- c("dataset", "K", "CV_error", "LogLikelihood")

label_map <- c(
  rice_pruned     = "S0: No mind (413)",   rice_pruned2    = "S1: mind 0.10 after QC (409)",
  rice_pruned3    = "S3: mind 0.10 before QC (379)", sens_D1_S2 = "S2: mind 0.05 (293)",
  sens_D1_S3      = "S4: mind 0.02 (98)",
  sens_D2_M0      = "M0: MAF none (1638)", sens_D2_M1  = "M1: MAF 0.01 (1601)",
  sens_D2_M2      = "M2: MAF 0.03 (1386)", sens_D2_M4  = "M4: MAF 0.10 (844)",
  sens_D3_Gnone   = "G0: GENO none (1437)", sens_D3_G010 = "G2: GENO 0.10 (1299)",
  sens_D3_G002    = "G4: GENO 0.02 (982)",
  sens_D4_LDnone  = "L0: No pruning (26474)", sens_D4_LD08  = "L1: r2<0.8 (9609)",
  sens_D4_LD05    = "L2: r2<0.5 (4278)",   sens_D4_LD100 = "L4: 100/10/0.2 (822)"
)

# ---- FIGURE 1: WORKFLOW PIPELINE ----
pdf("Figure_5_workflow.pdf", width=10, height=8)
par(mar=c(0,0,0,0))
plot(0, type="n", xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")

# Define boxes
boxes <- list(
  list(label="Raw data\n413 accessions\n36,901 SNPs",     x1=0.05, y1=0.85, x2=0.28, y2=0.97, col="#FBB4AE"),
  list(label="Sample QC\n--mind 0.10\n383 accessions",     x1=0.38, y1=0.85, x2=0.61, y2=0.97, col="#B3CDE3"),
  list(label="Remove 4 outliers\n379 accessions",          x1=0.71, y1=0.85, x2=0.94, y2=0.97, col="#CCEBC5"),
  list(label="Marker QC\n--geno 0.05\n26,474 SNPs",        x1=0.38, y1=0.65, x2=0.61, y2=0.77, col="#B3CDE3"),
  list(label="MAF filter\n--maf 0.05\n26,474 SNPs",        x1=0.38, y1=0.45, x2=0.61, y2=0.57, col="#B3CDE3"),
  list(label="LD pruning\n--indep-pairwise 50 5 0.2\n1,187 SNPs", x1=0.38, y1=0.25, x2=0.61, y2=0.37, col="#DECBE4"),
  list(label="Final dataset\n379 x 1,187",                 x1=0.05, y1=0.05, x2=0.28, y2=0.17, col="#FED9A6"),
  list(label="PCA\nPopulation\nstructure",                 x1=0.38, y1=0.05, x2=0.52, y2=0.17, col="#FFFFCC"),
  list(label="ADMIXTURE\nK=1-12\nAncestry\nproportions",   x1=0.56, y1=0.05, x2=0.70, y2=0.17, col="#FFFFCC"),
  list(label="Diversity\nObs/Exp Het\nFIS",                x1=0.74, y1=0.05, x2=0.88, y2=0.17, col="#FFFFCC")
)

# Draw boxes
for (b in boxes) {
  rect(b$x1, b$y1, b$x2, b$y2, col=b$col, border="gray40", lwd=1.5, xpd=NA)
  text((b$x1+b$x2)/2, (b$y1+b$y2)/2, b$label, cex=0.7, font=2, xpd=NA)
}

# Arrows
arrows(0.28, 0.91, 0.37, 0.91, length=0.08, lwd=2, xpd=NA)
arrows(0.61, 0.91, 0.70, 0.91, length=0.08, lwd=2, xpd=NA)
arrows(0.495, 0.77, 0.495, 0.68, length=0.08, lwd=2, xpd=NA)
arrows(0.495, 0.57, 0.495, 0.48, length=0.08, lwd=2, xpd=NA)
arrows(0.495, 0.37, 0.495, 0.28, length=0.08, lwd=2, xpd=NA)
arrows(0.28, 0.11, 0.37, 0.11, length=0.08, lwd=2, xpd=NA)

text(0.18, 0.78, "QC pipeline", cex=0.8, font=3, col="gray40", xpd=NA)
text(0.18, 0.40, "Final", cex=0.8, font=3, col="gray40", xpd=NA)
text(0.63, 0.40, "Downstream analyses", cex=0.8, font=3, col="gray40", xpd=NA)
dev.off()
cat("Figure 5 (workflow) saved\n")

# ---- FIGURE 2: QC IMPACT ----
master <- read.csv("Table_S4_sensitivity_master.csv", stringsAsFactors=FALSE)

# Add baseline (rice_pruned3, the reference) row
baseline_row <- master[master$ID == "rice_pruned3", ]

# Organize by dimension for barplots
dim_order <- list(
  "Sample QC" = c("S0: No mind (413)", "S1: mind 0.10 after QC (409)",
                  "S3: mind 0.10 before QC (379)", "S2: mind 0.05 (293)", "S4: mind 0.02 (98)"),
  "MAF" = c("M0: MAF none (1638)", "M1: MAF 0.01 (1601)", "M2: MAF 0.03 (1386)",
            "M4: MAF 0.10 (844)"),
  "GENO" = c("G0: GENO none (1437)", "G2: GENO 0.10 (1299)", "G4: GENO 0.02 (982)"),
  "LD pruning" = c("L0: No pruning (26474)", "L1: r2<0.8 (9609)",
                   "L2: r2<0.5 (4278)", "L4: 100/10/0.2 (822)")
)

# Metrics to plot, with normalization to baseline
metrics <- list(
  "Samples" = list(col="steelblue", digits=0),
  "SNPs" = list(col="forestgreen", digits=0),
  "CV_error" = list(col="tomato", digits=3),
  "Fst" = list(col="purple", digits=3),
  "PC1_corr" = list(col="darkorange", digits=3)
)

pdf("Figure_6_QC_impact.pdf", width=14, height=10)
par(mfrow=c(2,3), mar=c(6,4,3,1), mgp=c(2.5,0.7,0))

# Plot 1: Samples retained
plot(NULL, xlim=c(0.5,5.5), ylim=c(0,450), xaxt="n",
     xlab="", ylab="Samples retained", main="A) Sample filtering impact", cex.main=1)
axis(1, at=1:5, labels=dim_order[["Sample QC"]], las=2, cex.axis=0.55)
bl <- 379
vals <- c(413, 409, bl, 293, 98)
cols <- ifelse(vals < bl, "tomato", ifelse(vals > bl, "forestgreen", "gray40"))
barplot(vals, add=TRUE, axes=FALSE, col=cols, names.arg=rep("",5))
abline(h=bl, lty=2, col="gray40", lwd=1.5)
text(0.7, bl+10, paste0("baseline\n(", bl, ")"), cex=0.5, col="gray40")

# Plot 2: SNPs retained (MAF dimension)
plot(NULL, xlim=c(0.5,4.5), ylim=c(0,3000), xaxt="n",
     xlab="", ylab="SNPs retained", main="B) MAF threshold impact", cex.main=1)
axis(1, at=1:4, labels=dim_order[["MAF"]], las=2, cex.axis=0.55)
bl <- 1187
vals <- c(1638, 1601, 1386, 844)
cols <- ifelse(vals < bl, "tomato", ifelse(vals > bl, "forestgreen", "gray40"))
barplot(vals, add=TRUE, axes=FALSE, col=cols, names.arg=rep("",4))
abline(h=bl, lty=2, col="gray40", lwd=1.5)
text(0.7, bl+50, paste0("baseline\n(", bl, ")"), cex=0.5, col="gray40")

# Plot 3: CV error at K=5
plot(NULL, xlim=c(0.5,4.5), ylim=c(0.4,0.8), xaxt="n",
     xlab="", ylab="CV error (K=5)", main="C) CV error across dimensions", cex.main=1)
axis(1, at=1:4, labels=names(dim_order), las=1, cex.axis=0.7)
bl_cv <- master$CV_K5[master$ID=="rice_pruned3"]
for (i in seq_along(dim_order)) {
  dn <- names(dim_order)[i]
  labels <- dim_order[[i]]
  vals <- sapply(labels, function(l) {
    master$CV_K5[master$Label %in% l]
  })
  points(rep(i, length(vals)), vals, pch=19, col=i, cex=1.2)
  segments(i, min(vals), i, max(vals), col=i, lwd=2)
}
abline(h=bl_cv, lty=2, col="gray40", lwd=1.5)
text(4.3, bl_cv+0.005, paste0("baseline\nCV=", round(bl_cv,3)), cex=0.5, col="gray40")
legend("bottomleft", legend=names(dim_order), pch=19, col=1:4, cex=0.6)

# Plot 4: FST
plot(NULL, xlim=c(0.5,4.5), ylim=c(0.2,0.65), xaxt="n",
     xlab="", ylab="Mean Fst", main="D) Fst across dimensions", cex.main=1)
axis(1, at=1:4, labels=names(dim_order), las=1, cex.axis=0.7)
bl_fst <- master$Fst[master$ID=="rice_pruned3"]
for (i in seq_along(dim_order)) {
  dn <- names(dim_order)[i]
  labels <- dim_order[[i]]
  vals <- sapply(labels, function(l) {
    master$Fst[master$Label %in% l]
  })
  points(rep(i, length(vals)), vals, pch=19, col=i, cex=1.2)
  segments(i, min(vals), i, max(vals), col=i, lwd=2)
}
abline(h=bl_fst, lty=2, col="gray40", lwd=1.5)
text(4.3, bl_fst+0.005, paste0("baseline\nFst=", round(bl_fst,3)), cex=0.5, col="gray40")

# Plot 5: PC1 correlation
master_pc <- master[!is.na(master$PC1_corr), ]
plot(NULL, xlim=c(0.5, nrow(master_pc)+0.5), ylim=c(0.85,1), xaxt="n",
     xlab="", ylab="PC1 correlation with baseline", main="E) PCA concordance", cex.main=1)
axis(1, at=1:nrow(master_pc), labels=master_pc$Label, las=2, cex.axis=0.45)
barplot(master_pc$PC1_corr, add=TRUE, axes=FALSE, col="darkorange",
        names.arg=rep("", nrow(master_pc)))
abline(h=0.99, lty=2, col="red", lwd=1.5)
text(nrow(master_pc)-0.5, 0.991, "r = 0.99", cex=0.6, col="red")

# Plot 6: Legend / summary
plot.new()
legend("center", legend=c("Above baseline (more conservative)",
                          "Below baseline (more liberal)",
                          "Baseline value"), 
       fill=c("forestgreen", "tomato", "gray40"), cex=0.9, bty="n")
text(0.5, 0.3, "Baseline: mind 0.10 + MAF 0.05 +\nGENO 0.05 + LD r2<0.2\n379 samples, 1187 SNPs", cex=0.8)
dev.off()
cat("Figure 6 (QC impact) saved\n")

# ---- FIGURE 3: CLUSTER STABILITY ----
scenarios_k5 <- list(
  Sample_S0 = "rice_pruned", Sample_S1 = "rice_pruned2",
  Sample_S3 = "rice_pruned3", Sample_S2 = "sens_D1_S2", Sample_S4 = "sens_D1_S3",
  MAF_M0 = "sens_D2_M0", MAF_M1 = "sens_D2_M1", MAF_M2 = "sens_D2_M2", MAF_M4 = "sens_D2_M4",
  Geno_G0 = "sens_D3_Gnone", Geno_G2 = "sens_D3_G010", Geno_G4 = "sens_D3_G002",
  LD_L0 = "sens_D4_LDnone", LD_L1 = "sens_D4_LD08", LD_L2 = "sens_D4_LD05", LD_L4 = "sens_D4_LD100"
)

# Read Q files and assign clusters
read_q <- function(ds) {
  qfile <- paste0(SENS, "/", ds, "_K5.Q")
  if (!file.exists(qfile)) return(NULL)
  q <- as.matrix(read.table(qfile))
  apply(q, 1, which.max)
}

cat("Computing cluster stability...\n")
all_clusters <- list()
fam_list <- list()

for (nm in names(scenarios_k5)) {
  ds <- scenarios_k5[[nm]]
  cl <- read_q(ds)
  if (is.null(cl)) next
  all_clusters[[nm]] <- cl
  
  fam <- read.table(paste0(BD, "/", ds, ".fam"))
  fam_list[[nm]] <- fam
}

# Find samples common to all scenarios
# Start with samples from the baseline (rice_pruned3)
baseline_sc <- names(scenarios_k5)[grep("Sample_S3", names(scenarios_k5))]
bl_samples <- fam_list[[baseline_sc]][,2]
bl_fid <- fam_list[[baseline_sc]][,1]
bl_keys <- paste0(bl_fid, "_", bl_samples)

# Build a matrix: samples x scenarios
n_scenarios <- length(all_clusters)
sc_names <- names(all_clusters)

# Use baseline sample set and find those present in each scenario
stability <- data.frame(FID=bl_fid, IID=bl_samples, stringsAsFactors=FALSE)

for (nm in sc_names) {
  fam_df <- fam_list[[nm]]
  keys <- paste0(fam_df[,1], "_", fam_df[,2])
  
  # Match baseline samples to this scenario
  clust_col <- rep(NA, nrow(stability))
  names(clust_col) <- bl_keys
  match_idx <- match(bl_keys, keys)
  valid <- !is.na(match_idx)
  clust_col[valid] <- all_clusters[[nm]][match_idx[valid]]
  stability[[nm]] <- clust_col
}

# Align cluster labels across scenarios (exhaustive search, K <= 6)
ref_col <- grep("Sample_S3", sc_names, value=TRUE)
ref <- stability[[ref_col]]

align_labels <- function(ref_vec, test_vec) {
  valid <- !is.na(ref_vec) & !is.na(test_vec)
  ref_v <- ref_vec[valid]
  test_v <- test_vec[valid]
  
  n_clust <- max(c(ref_v, test_v), na.rm=TRUE)
  if (n_clust > 6) return(1:n_clust)
  
  # Build confusion matrix: rows = test clusters, cols = ref clusters
  cm <- matrix(0, n_clust, n_clust)
  for (i in seq_along(ref_v)) {
    cm[test_v[i], ref_v[i]] <- cm[test_v[i], ref_v[i]] + 1
  }
  
  # Greedy matching: assign each test cluster to its best matching ref cluster
  used_ref <- rep(FALSE, n_clust)
  mapping <- rep(NA, n_clust)
  for (iter in 1:n_clust) {
    best_val <- -1
    best_t <- 0
    best_r <- 0
    for (t in 1:n_clust) {
      if (!is.na(mapping[t])) next
      for (r in 1:n_clust) {
        if (used_ref[r]) next
        if (cm[t, r] > best_val) {
          best_val <- cm[t, r]
          best_t <- t
          best_r <- r
        }
      }
    }
    mapping[best_t] <- best_r
    used_ref[best_r] <- TRUE
  }
  mapping
}

# Apply label alignment to all non-reference scenarios
for (nm in sc_names) {
  if (nm == ref_col) next
  test_vec <- stability[[nm]]
  mapping <- align_labels(ref, test_vec)
  
  # Apply mapping
  aligned <- test_vec
  valid <- !is.na(test_vec)
  aligned[valid] <- mapping[test_vec[valid]]
  stability[[nm]] <- aligned
  
  # Verify alignment quality
  valid_both <- !is.na(ref) & !is.na(aligned)
  agree_rate <- mean(ref[valid_both] == aligned[valid_both], na.rm=TRUE)
  cat("  Alignment", nm, ": agreement rate =", round(agree_rate, 3), "\n")
}

# Now compute per-sample agreement with baseline
ref <- stability[[ref_col]]
agreement <- sapply(1:nrow(stability), function(i) {
  if (is.na(ref[i])) return(NA)
  other_sc <- sc_names[sc_names != ref_col]
  others <- unlist(stability[i, other_sc])
  sum(others == ref[i], na.rm=TRUE) / sum(!is.na(others))
})

stability$Agreement <- agreement

# Distribution of stability
cat("Cluster stability summary:\n")
cat("  Mean agreement:", round(mean(agreement, na.rm=TRUE), 3), "\n")
cat("  Median agreement:", round(median(agreement, na.rm=TRUE), 3), "\n")
cat("  % fully stable (agree=1):", round(mean(agreement >= 0.999, na.rm=TRUE)*100, 1), "%\n")
cat("  % highly stable (agree>=0.9):", round(mean(agreement >= 0.9, na.rm=TRUE)*100, 1), "%\n")
cat("  % unstable (agree<0.8):", round(mean(agreement < 0.8, na.rm=TRUE)*100, 1), "%\n")

# Count how many scenarios each sample appears in
stability$N_scenarios <- rowSums(!is.na(stability[, sc_names]))

# Write stability table
write.csv(stability, "Table_S5_cluster_stability.csv", row.names=FALSE)
cat("Table S5 (cluster stability) saved\n")

# Stability barplot
pdf("Figure_7_cluster_stability.pdf", width=12, height=6)
par(mfrow=c(1,2), mar=c(4,4,3,1))

# Histogram of stability
valid_s <- stability[!is.na(stability$Agreement) & stability$N_scenarios >= 10, ]
hist(valid_s$Agreement, breaks=20, col="steelblue", border="white",
     xlab="Proportion of scenarios with same assignment", 
     ylab="Number of accessions",
     main=paste0("Cluster stability across filtering scenarios\n(n=", nrow(valid_s), " accessions, ", n_scenarios, " scenarios)"),
     cex.main=0.9)
abline(v=mean(valid_s$Agreement, na.rm=TRUE), lty=2, col="red", lwd=2)
text(mean(valid_s$Agreement, na.rm=TRUE)+0.02, par("usr")[4]*0.9,
     paste0("Mean=", round(mean(valid_s$Agreement, na.rm=TRUE), 3)),
     col="red", cex=0.7)

# Barplot of agreement by cluster (using baseline cluster)
bl_clusters <- stability[[ref_col]]
agree_by_cluster <- tapply(stability$Agreement, bl_clusters, function(x) {
  c(mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE), n=sum(!is.na(x)))
})
agree_mat <- do.call(rbind, agree_by_cluster)

bp <- barplot(agree_mat[, "mean"], ylim=c(0,1),
              col=rainbow(nrow(agree_mat)),
              xlab="Baseline cluster (K=5)", ylab="Mean agreement",
              main="Stability by baseline cluster assignment",
              cex.main=0.9, names.arg=paste0("K", 1:nrow(agree_mat)))
arrows(bp, agree_mat[, "mean"] - agree_mat[, "sd"],
       bp, agree_mat[, "mean"] + agree_mat[, "sd"],
       angle=90, code=3, length=0.05)
text(bp, agree_mat[, "mean"] + 0.03, paste0("n=", agree_mat[, "n"]), cex=0.6)
dev.off()
cat("Figure 7 (cluster stability) saved\n")

cat("\n===== ALL FIGURES COMPLETE =====\n")
