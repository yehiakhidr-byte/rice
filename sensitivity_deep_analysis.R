# Deep sensitivity analysis: cluster, PCA, diversity, FST, AMOVA, LD decay

library(adegenet)
library(vegan)
library(ape)
library(poppr)
library(cluster)

SENS_DIR <- "/home/yehia/sensitivity"
OUT_DIR  <- "/mnt/d/project genomic analysis"

# ---- 0. Define scenarios ----
scenarios <- c(
  "sens_D1_S2", "sens_D1_S3",
  "sens_D2_M0", "sens_D2_M1", "sens_D2_M2", "sens_D2_M4",
  "sens_D3_G002", "sens_D3_G010", "sens_D3_G020", "sens_D3_Gnone",
  "sens_D4_LD05", "sens_D4_LD08", "sens_D4_LD100", "sens_D4_LDnone",
  "rice_pruned", "rice_pruned2", "rice_pruned3"
)

pretty_names <- c(
  "S2: mind=0.05", "S3: mind=0.10",
  "M0: MAF=none", "M1: MAF=0.01", "M2: MAF=0.03", "M4: MAF=0.10",
  "G2: GENO=0.02", "G1: GENO=0.10", "G0: GENO=0.20", "Gnone: GENO=none",
  "L_0.5: LD<0.5", "L_0.8: LD<0.8", "L100: 100/10/0.2", "Lnone: no LD",
  "rice_pruned", "rice_pruned2", "rice_pruned3"
)
names(pretty_names) <- scenarios

cluster_names <- c("indica", "japonica", "aus", "aromatic", "admixed")

# ---- 1. Read Q files and FAM files ----
read_q <- function(file, ncols=5) {
  if (!file.exists(file)) return(NULL)
  m <- as.matrix(read.table(file, header=FALSE))
  max_idx <- apply(m, 1, which.max)
  attr(max_idx, "Q") <- m
  max_idx
}

read_clusters <- function(sc) {
  qfile <- file.path(SENS_DIR, paste0(sc, "_K5.Q"))
  famfile <- file.path(SENS_DIR, paste0(sc, ".fam"))
  if (!file.exists(qfile) || !file.exists(famfile)) return(NULL)
  fam <- read.table(famfile, header=FALSE)
  q <- read_q(qfile)
  list(
    sc = sc,
    samples = as.character(fam[,2]),
    cluster = q,
    Q = attr(q, "Q"),
    n = length(q),
    label = cluster_names[q]
  )
}

all_clusters <- lapply(scenarios, read_clusters)
names(all_clusters) <- scenarios
all_clusters <- all_clusters[!sapply(all_clusters, is.null)]

cat("Loaded", length(all_clusters), "scenario clusterings\n\n")

# ---- Identify baseline ----
baseline_sc <- "rice_pruned3"
cat("Baseline scenario:", baseline_sc, "\n")

# ---- 2. Confusion matrices per dimension ----
cat("\n===== 1. CONFUSION MATRICES =====\n")

# For each dimension, find matching scenarios and build confusion matrix vs baseline
dimensions <- list(
  Sample_filtering = c("sens_D1_S2", "sens_D1_S3"),
  MAF = c("sens_D2_M0", "sens_D2_M1", "sens_D2_M2", "sens_D2_M4"),
  GENO = c("sens_D3_G002", "sens_D3_G010", "sens_D3_G020", "sens_D3_Gnone"),
  LD_pruning = c("sens_D4_LD05", "sens_D4_LD08", "sens_D4_LD100", "sens_D4_LDnone")
)

build_confusion <- function(sc1, sc2, name1, name2) {
  d1 <- all_clusters[[sc1]]; d2 <- all_clusters[[sc2]]
  if (is.null(d1) || is.null(d2)) return(NULL)
  common <- intersect(d1$samples, d2$samples)
  idx1 <- match(common, d1$samples)
  idx2 <- match(common, d2$samples)
  c1 <- d1$cluster[idx1]; c2 <- d2$cluster[idx2]
  tab <- table(Baseline=cluster_names[c1], New=cluster_names[c2])
  attr(tab, "n_common") <- length(common)
  attr(tab, "n_mismatch") <- sum(c1 != c2)
  tab
}

confusion_list <- list()

for (dim_name in names(dimensions)) {
  sc_list <- dimensions[[dim_name]]
  cat("\n---", dim_name, "---\n")
  for (sc in sc_list) {
    if (all_clusters[[sc]]$n < 10) { cat(sc, "- too few samples\n"); next }
    cm <- build_confusion(baseline_sc, sc, "Baseline", sc)
    if (!is.null(cm)) {
      mismatch_pct <- round(100 * (attr(cm, "n_common") - sum(diag(cm))) / attr(cm, "n_common"), 1)
      cat("  vs", pretty_names[sc], ": n=", attr(cm, "n_common"), 
          ", mismatches=", mismatch_pct, "%\n", sep="")
      print(cm)
      confusion_list[[paste(baseline_sc, sc, sep="_vs_")]] <- cm
    }
  }
}

# ---- 3. Pairwise ARI / NMI ----
cat("\n\n===== 2. PAIRWISE ARI/NMI =====\n")

# Implement ARI
calculate_ari <- function(c1, c2) {
  valid <- c1 > 0 & c2 > 0
  c1 <- c1[valid]; c2 <- c2[valid]
  if (length(c1) < 5) return(NA)
  
  tab <- table(c1, c2)
  n <- sum(tab)
  if (n < 5) return(NA)
  
  sum_comb_ij <- sum(choose(tab, 2))
  sum_comb_i <- sum(choose(rowSums(tab), 2))
  sum_comb_j <- sum(choose(colSums(tab), 2))
  comb_n <- choose(n, 2)
  
  expected <- sum_comb_i * sum_comb_j / comb_n
  max_val <- (sum_comb_i + sum_comb_j) / 2
  
  (sum_comb_ij - expected) / (max_val - expected)
}

# Implement NMI
calculate_nmi <- function(c1, c2) {
  valid <- c1 > 0 & c2 > 0
  c1 <- c1[valid]; c2 <- c2[valid]
  if (length(c1) < 5) return(NA)
  
  tab <- table(c1, c2)
  n <- sum(tab)
  if (n < 5) return(NA)
  
  # Mutual information
  mi <- 0
  for (i in 1:nrow(tab)) {
    for (j in 1:ncol(tab)) {
      if (tab[i,j] > 0) {
        mi <- mi + (tab[i,j]/n) * log((tab[i,j] * n) / (rowSums(tab)[i] * colSums(tab)[j]))
      }
    }
  }
  
  # Entropies
  h1 <- -sum((rowSums(tab)/n) * log(rowSums(tab)/n))
  h2 <- -sum((colSums(tab)/n) * log(colSums(tab)/n))
  
  if (h1 == 0 || h2 == 0) return(NA)
  mi / sqrt(h1 * h2)
}

# Find scenarios that share samples for pairwise comparison
# All sens_D* datasets have 383 samples; rice_* may differ
sens_names <- names(all_clusters)
sens_names <- sens_names[sapply(all_clusters, function(x) x$n >= 50)]

n_sc <- length(sens_names)
ari_matrix <- matrix(NA, n_sc, n_sc, dimnames=list(pretty_names[sens_names], pretty_names[sens_names]))
nmi_matrix <- matrix(NA, n_sc, n_sc, dimnames=list(pretty_names[sens_names], pretty_names[sens_names]))

for (i in 1:n_sc) {
  for (j in i:n_sc) {
    s1 <- sens_names[i]; s2 <- sens_names[j]
    d1 <- all_clusters[[s1]]; d2 <- all_clusters[[s2]]
    common <- intersect(d1$samples, d2$samples)
    if (length(common) < 10) next
    idx1 <- match(common, d1$samples)
    idx2 <- match(common, d2$samples)
    ari <- calculate_ari(d1$cluster[idx1], d2$cluster[idx2])
    nmi <- calculate_nmi(d1$cluster[idx1], d2$cluster[idx2])
    ari_matrix[i,j] <- ari_matrix[j,i] <- ari
    nmi_matrix[i,j] <- nmi_matrix[j,i] <- nmi
  }
}

# Find the baseline row for ARI
baseline_idx <- which(sens_names == baseline_sc)
if (length(baseline_idx) > 0) {
  cat("\nARI vs baseline (", pretty_names[baseline_sc], "):\n", sep="")
  ari_vs_baseline <- ari_matrix[baseline_idx, ]
  print(round(ari_vs_baseline[!is.na(ari_vs_baseline)], 3))
}

cat("\nNMI vs baseline:\n")
nmi_vs_baseline <- nmi_matrix[baseline_idx, ]
print(round(nmi_vs_baseline[!is.na(nmi_vs_baseline)], 3))

# Write ARI/NMI tables
write.csv(ari_matrix, file.path(OUT_DIR, "Table_S6_ARI_matrix.csv"))
write.csv(nmi_matrix, file.path(OUT_DIR, "Table_S7_NMI_matrix.csv"))

# Heatmap for ARI
pdf(file.path(OUT_DIR, "Figure_8_ARI_heatmap.pdf"), width=10, height=9)
par(mar=c(10,10,3,2))
heatmap(ari_matrix, symm=TRUE, col=colorRampPalette(c("white","yellow","orange","red"))(100),
        main="Adjusted Rand Index between scenarios", cexRow=0.6, cexCol=0.6)
dev.off()
cat("Figure 8 (ARI heatmap) saved\n")

pdf(file.path(OUT_DIR, "Figure_9_NMI_heatmap.pdf"), width=10, height=9)
par(mar=c(10,10,3,2))
heatmap(nmi_matrix, symm=TRUE, col=colorRampPalette(c("white","lightblue","blue","darkblue"))(100),
        main="Normalized Mutual Information between scenarios", cexRow=0.6, cexCol=0.6)
dev.off()
cat("Figure 9 (NMI heatmap) saved\n")

# ---- 4. PCA correlation ----
cat("\n\n===== 3. PCA SIMILARITY =====\n")

read_pca <- function(sc) {
  f <- file.path(SENS_DIR, paste0(sc, "_pca.eigenvec"))
  if (!file.exists(f)) return(NULL)
  d <- read.table(f, header=FALSE)
  list(samples=as.character(d[,2]), PC1=d[,3], PC2=d[,4])
}

pca_list <- lapply(sens_names, read_pca)
names(pca_list) <- sens_names
pca_list <- pca_list[!sapply(pca_list, is.null)]

# PC correlation vs baseline
baseline_pca <- pca_list[[baseline_sc]]
if (!is.null(baseline_pca)) {
  cat("PC correlation vs baseline (", pretty_names[baseline_sc], "):\n", sep="")
  pc_cor <- data.frame(Scenario=character(), PC1_r=numeric(), PC2_r=numeric(), stringsAsFactors=FALSE)
  for (sc in names(pca_list)) {
    if (sc == baseline_sc) next
    dp <- pca_list[[sc]]
    if (is.null(dp)) next
    common <- intersect(baseline_pca$samples, dp$samples)
    if (length(common) < 10) next
    i1 <- match(common, baseline_pca$samples)
    i2 <- match(common, dp$samples)
    r1 <- cor(baseline_pca$PC1[i1], dp$PC1[i2])
    r2 <- cor(baseline_pca$PC2[i1], dp$PC2[i2])
    pc_cor <- rbind(pc_cor, data.frame(Scenario=pretty_names[sc], PC1_r=round(r1,3), PC2_r=round(r2,3)))
    cat("  ", pretty_names[sc], ": PC1 r=", round(r1,3), ", PC2 r=", round(r2,3), "\n", sep="")
  }
  write.csv(pc_cor, file.path(OUT_DIR, "Table_S8_PCA_correlation.csv"), row.names=FALSE)
}

# ---- 5. Genetic diversity ----
cat("\n\n===== 4. GENETIC DIVERSITY =====\n")

# PIC = 1 - sum(p_i^2) for biallelic marker = He
# Shannon = -sum(p_i * ln(p_i))
# Ho from HWE output
# He from freq output

read_freq <- function(sc) {
  f <- file.path(SENS_DIR, paste0(sc, "_freq.frq"))
  if (!file.exists(f)) return(NULL)
  frq <- read.table(f, header=TRUE, stringsAsFactors=FALSE)
  frq
}

read_hwe <- function(sc) {
  f <- file.path(SENS_DIR, paste0(sc, "_hwe.hwe"))
  if (!file.exists(f)) return(NULL)
  hwe <- read.table(f, header=TRUE, stringsAsFactors=FALSE)
  hwe
}

diversity_table <- data.frame(Scenario=character(), N_SNPs=integer(), N_Samples=integer(),
                              PIC=numeric(), He=numeric(), Ho=numeric(), Shannon=numeric(),
                              stringsAsFactors=FALSE)

for (sc in sens_names) {
  frq <- read_freq(sc)
  if (is.null(frq)) next
  n_snps <- nrow(frq)
  n_samp <- all_clusters[[sc]]$n
  
  # MAF = minor allele frequency
  maf <- ifelse(frq$MAF > 0.5, 1 - frq$MAF, frq$MAF)
  maj_f <- 1 - maf
  
  # PIC = 1 - sum(p_i^2) = 2*MAF*(1-MAF) for biallelic
  pic <- 2 * maf * maj_f  # same as He
  
  # Shannon Index = -sum(p_i * ln(p_i))
  shannon <- -(maf * log(maf + 1e-10) + maj_f * log(maj_f + 1e-10))
  
  mean_pic <- mean(pic, na.rm=TRUE)
  mean_he <- mean(pic, na.rm=TRUE)  # He = PIC for biallelic
  mean_shannon <- mean(shannon, na.rm=TRUE)
  
  # Ho from HWE (O.HET. is count, need to compute rate from GENO column)
  hwe <- read_hwe(sc)
  if (!is.null(hwe)) {
    # Parse GENO column "N11/N12/N22" to get total N
    geno_parts <- strsplit(as.character(hwe$GENO), "/")
    n_total <- sapply(geno_parts, function(x) sum(as.numeric(x), na.rm=TRUE))
    ho_rate <- hwe$O.HET. / n_total
    mean_ho <- mean(ho_rate, na.rm=TRUE)
  } else {
    mean_ho <- NA
  }
  
  diversity_table <- rbind(diversity_table, data.frame(
    Scenario = pretty_names[sc], N_SNPs = n_snps, N_Samples = n_samp,
    PIC = round(mean_pic, 4), He = round(mean_he, 4), Ho = round(mean_ho, 4),
    Shannon = round(mean_shannon, 4), stringsAsFactors=FALSE
  ))
  cat("  ", pretty_names[sc], ": PIC=", round(mean_pic,4), 
      " He=", round(mean_he,4), " Ho=", round(mean_ho,4), 
      " Shannon=", round(mean_shannon,4), "\n", sep="")
}

write.csv(diversity_table, file.path(OUT_DIR, "Table_S9_genetic_diversity.csv"), row.names=FALSE)

# Diversity barplot
pdf(file.path(OUT_DIR, "Figure_10_diversity.pdf"), width=12, height=8)
par(mfrow=c(2,2), mar=c(8,4,3,1))
for (metric in c("PIC", "He", "Ho", "Shannon")) {
  cols <- rainbow(nrow(diversity_table))
  barplot(diversity_table[[metric]], names.arg=diversity_table$Scenario, las=2,
          col=cols, main=paste("Mean", metric), ylab=metric, cex.names=0.6)
  abline(h=mean(diversity_table[[metric]], na.rm=TRUE), lty=2, col="red")
}
dev.off()
cat("Figure 10 (diversity) saved\n")

# ---- 6. FST comparison ----
cat("\n\n===== 5. FST COMPARISON =====\n")

# Read existing FST summary
fst_file <- file.path(SENS_DIR, "fst_summary.csv")
if (file.exists(fst_file)) {
  fst <- read.csv(fst_file, stringsAsFactors=FALSE)
  cat("FST summary loaded:", nrow(fst), "scenarios\n")
  print(fst[, c("Scenario", "Samples", "SNPs", "Fst_mean", "Fst_sd")])
  
  # FST barplot
  pdf(file.path(OUT_DIR, "Figure_11_FST_comparison.pdf"), width=10, height=6)
  par(mar=c(8,4,3,1))
  fst_mean <- fst$Fst_mean
  fst_sd <- ifelse(is.na(fst$Fst_sd) | fst$Fst_sd == "", 0, as.numeric(fst$Fst_sd))
  bp <- barplot(fst_mean, names.arg=fst$Scenario, las=2, col="steelblue",
                main="Mean Pairwise FST across scenarios", ylab="Mean FST", cex.names=0.7,
                ylim=c(0, max(fst_mean + fst_sd, na.rm=TRUE)))
  arrows(bp, fst_mean - fst_sd, bp, fst_mean + fst_sd, angle=90, code=3, length=0.05)
  abline(h=mean(fst_mean, na.rm=TRUE), lty=2, col="red")
  dev.off()
  cat("Figure 11 (FST comparison) saved\n")
} else {
  cat("FST summary not found at", fst_file, "\n")
}

# ---- 7. AMOVA via PLINK distance + vegan::adonis2 ----
cat("\n\n===== 6. AMOVA =====\n")

base_clust <- all_clusters[[baseline_sc]]
if (!is.null(base_clust) && base_clust$n >= 10) {
  base_samples <- base_clust$samples
  base_pops <- base_clust$cluster

  amova_results <- data.frame(Scenario=character(), N_Samples=integer(),
                              R2=numeric(), F_Stat=numeric(), Phi_ST=numeric(),
                              stringsAsFactors=FALSE)

  for (sc in sens_names) {
    scd <- all_clusters[[sc]]
    if (is.null(scd) || scd$n < 10) next

    common <- intersect(scd$samples, base_samples)
    if (length(common) < 20) next
    base_idx <- match(common, base_samples)
    sc_idx   <- match(common, scd$samples)
    pop_labels <- cluster_names[base_pops[base_idx]]

    mdist_file <- file.path(SENS_DIR, paste0(sc, "_dist.mdist"))
    mdist_id_file <- file.path(SENS_DIR, paste0(sc, "_dist.mdist.id"))
    if (!file.exists(mdist_file) || !file.exists(mdist_id_file)) next

    ids <- read.table(mdist_id_file, header=FALSE, stringsAsFactors=FALSE)[,2]
    id_common <- intersect(ids, common)
    if (length(id_common) < 20) next

    m <- as.matrix(read.table(mdist_file, header=FALSE))
    idx <- match(id_common, ids)
    m_sub <- m[idx, idx]
    colnames(m_sub) <- rownames(m_sub) <- id_common

    pop_f <- factor(pop_labels[match(id_common, common)])

    set.seed(999)
    ad <- adonis2(as.dist(m_sub) ~ pop_f, permutations=999)
    # adonis2 row name is "Model" for simple formulas
    r2 <- ad[1, "R2"]
    f_stat <- ad[1, "F"]
    phi_st <- r2  # For 1-level AMOVA with Euclidean distance, R2 = Phi_ST

    amova_results <- rbind(amova_results, data.frame(
      Scenario = pretty_names[sc],
      N_Samples = length(id_common),
      R2 = round(r2, 4),
      F_Stat = round(f_stat, 2),
      Phi_ST = round(phi_st, 4),
      stringsAsFactors=FALSE
    ))
    cat("  ", pretty_names[sc], ": Phi_ST =", round(phi_st, 4), " F =", round(f_stat, 2), "\n")
  }

  if (nrow(amova_results) > 0) {
    write.csv(amova_results, file.path(OUT_DIR, "Table_S10_AMOVA.csv"), row.names=FALSE)
    cat("AMOVA results saved\n")
    print(amova_results)

    if (nrow(amova_results) > 1) {
      pdf(file.path(OUT_DIR, "Figure_12_AMOVA.pdf"), width=10, height=6)
      par(mar=c(8,4,3,1))
      cols <- ifelse(amova_results$Phi_ST >= mean(amova_results$Phi_ST, na.rm=TRUE), "darkgreen", "lightgreen")
      bp <- barplot(amova_results$Phi_ST, names.arg=amova_results$Scenario, las=2,
                    col=cols, main="Phi_ST (AMOVA - among-population variance)",
                    ylab="Phi_ST", cex.names=0.7)
      abline(h=mean(amova_results$Phi_ST, na.rm=TRUE), lty=2, col="red")
      dev.off()
      cat("Figure 12 (AMOVA) saved\n")
    }
  }
} else {
  cat("Baseline clustering not available for AMOVA\n")
}

# ---- 8. LD decay ----
cat("\n\n===== 7. LD DECAY =====\n")

ld_scenarios <- c("sens_D4_LDnone", "sens_D4_LD08", "sens_D4_LD05", "sens_D4_LD100")
ld_names <- c("No LD pruning", "r2 < 0.8", "r2 < 0.5", "Window 100/10, r2<0.2")
ld_colors <- c("black", "red", "blue", "darkgreen")

# Define distance bins (kb)
bins <- c(0, 10, 50, 100, 200, 500, 1000)
bin_labels <- c("0-10", "10-50", "50-100", "100-200", "200-500", "500-1000")

ld_decay_data <- list()

for (i in seq_along(ld_scenarios)) {
  sc <- ld_scenarios[i]
  ld_file <- file.path(SENS_DIR, paste0(sc, "_ld.ld"))
  if (!file.exists(ld_file)) {
    cat("  LD file not found:", ld_file, "\n")
    next
  }
  
  cat("  Reading", ld_names[i], "... ")
  ld <- read.table(ld_file, header=TRUE, stringsAsFactors=FALSE)
  cat(nrow(ld), "pairs\n")
  
  # Compute distance in kb from BP positions
  ld$DIST_KB <- abs(ld$BP_B - ld$BP_A) / 1000
  
  # Bin by distance
  ld$BIN <- cut(ld$DIST_KB, breaks=bins, labels=bin_labels, right=FALSE)
  
  # Mean r2 per bin
  bin_mean <- tapply(ld$R2, ld$BIN, mean, na.rm=TRUE)
  bin_sd <- tapply(ld$R2, ld$BIN, sd, na.rm=TRUE)
  bin_n <- tapply(ld$R2, ld$BIN, function(x) sum(!is.na(x)))
  
  # Count SNPs from BIM file
  bim_file <- file.path(SENS_DIR, paste0(sc, ".bim"))
  n_snps <- ifelse(file.exists(bim_file), as.integer(nrow(read.table(bim_file, header=FALSE))), NA)
  
  ld_decay_data[[sc]] <- list(
    name = ld_names[i],
    bins = bin_labels,
    mean_r2 = as.numeric(bin_mean),
    sd_r2 = as.numeric(bin_sd),
    n_pairs = as.numeric(bin_n),
    n_snps = n_snps,
    total_pairs = nrow(ld)
  )
  
  cat("    Bins:", paste(round(bin_mean, 3), collapse=", "), "\n")
}

# Plot LD decay curves
if (length(ld_decay_data) > 0) {
  pdf(file.path(OUT_DIR, "Figure_13_LD_decay.pdf"), width=8, height=6)
  par(mar=c(5,5,3,2))
  
  x_pos <- 1:length(bin_labels)
  
  plot(x_pos, type="n", xlim=c(0.5, length(bin_labels)+0.5), ylim=c(0, 0.5),
       xaxt="n", xlab="Physical distance (kb)", ylab="Mean r2",
       main="LD Decay by Pruning Scenario", cex.lab=1.2)
  axis(1, at=x_pos, labels=bin_labels)
  
  for (i in seq_along(ld_decay_data)) {
    dd <- ld_decay_data[[i]]
    if (length(dd$mean_r2) == 0) next
    lines(x_pos, dd$mean_r2, type="o", pch=16, col=ld_colors[i], lwd=2)
  }
  
  legend("topright", legend=names(ld_decay_data), col=ld_colors[1:length(ld_decay_data)], 
         lwd=2, pch=16, cex=0.8)
  dev.off()
  cat("Figure 13 (LD decay) saved\n")
  
  # Also create an LD decay data table
  ld_table <- data.frame(Bin=bin_labels, stringsAsFactors=FALSE)
  for (sc in names(ld_decay_data)) {
    dd <- ld_decay_data[[sc]]
    vals <- rep(NA, length(bin_labels))
    vals[1:length(dd$mean_r2)] <- round(dd$mean_r2, 4)
    ld_table[[dd$name]] <- vals
  }
  write.csv(ld_table, file.path(OUT_DIR, "Table_S11_LD_decay.csv"), row.names=FALSE)
  
  # Also save N_SNPs per scenario
  snp_info <- data.frame(Scenario=names(ld_decay_data), 
                         N_SNPs=sapply(ld_decay_data, `[[`, "n_snps"),
                         Total_Pairs=sapply(ld_decay_data, `[[`, "total_pairs"))
  write.csv(snp_info, file.path(OUT_DIR, "Table_S11b_LD_snp_counts.csv"), row.names=FALSE)
}

# ---- 9. Summary table ----
cat("\n\n===== SUMMARY =====\n")
cat("All deep sensitivity analyses complete.\n")
cat("Output files:\n")
cat("  Table_S6_ARI_matrix.csv\n")
cat("  Table_S7_NMI_matrix.csv\n")
cat("  Table_S8_PCA_correlation.csv\n")
cat("  Table_S9_genetic_diversity.csv\n")
cat("  Table_S10_AMOVA.csv\n")
cat("  Table_S11_LD_decay.csv\n")
cat("  Figure_8_ARI_heatmap.pdf\n")
cat("  Figure_9_NMI_heatmap.pdf\n")
cat("  Figure_10_diversity.pdf\n")
cat("  Figure_11_FST_comparison.pdf\n")
cat("  Figure_12_AMOVA.pdf\n")
cat("  Figure_13_LD_decay.pdf\n")
