library(adegenet)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(vegan)
library(poppr)
library(dplyr)

workdir <- "/home/yehia/sensitivity"
setwd(workdir)

# ============================================================
# 0. SCENARIO MAPPING
# ============================================================
scenarios <- data.frame(
  File = c("rice_pruned", "rice_pruned2", "sens_D1_S2", "rice_pruned3", "sens_D1_S3",
           "sens_D2_M0", "sens_D2_M1", "sens_D2_M2", "sens_D2_M4",
           "sens_D3_Gnone", "sens_D3_G020", "sens_D3_G010", "sens_D3_G002",
           "sens_D4_LDnone", "sens_D4_LD08", "sens_D4_LD05", "sens_D4_LD100"),
  Label = c("S0:NoMind", "S1:Mind10After", "S2:Mind05", "S3:Mind10Baseline", "S4:Mind02",
            "M0:MAFnone", "M1:MAF01", "M2:MAF03", "M4:MAF10",
            "G0:GENOnone", "G1:GENO20", "G2:GENO10", "G4:GENO02",
            "L0:LDnone", "L1:LD08", "L2:LD05", "L4:LD100"),
  Dim = c("Sample","Sample","Sample","Sample","Sample",
          "MAF","MAF","MAF","MAF",
          "GENO","GENO","GENO","GENO",
          "LD","LD","LD","LD"),
  Short = c("S0","S1","S2","S3","S4",
            "M0","M1","M2","M4",
            "G0","G1","G2","G4",
            "L0","L1","L2","L4"),
  stringsAsFactors = FALSE
)
rownames(scenarios) <- scenarios$Short
baseline_file <- "rice_pruned3"  # S3 = M3 = G3 = L3
baseline_label <- "S3:Mind10Baseline"

# ============================================================
# 1. READ AND ALIGN ALL CLUSTER ASSIGNMENTS
# ============================================================
cat("Reading Q files...\n")
read_q <- function(file, k=5) {
  q <- as.matrix(read.table(paste0(file, "_K5.Q"), header=FALSE))
  apply(q, 1, which.max)
}

assignments <- list()
for (i in 1:nrow(scenarios)) {
  f <- scenarios$File[i]
  lab <- scenarios$Label[i]
  assignments[[lab]] <- read_q(f)
}

# Baseline reference
ref_label <- "S3:Mind10Baseline"
ref <- assignments[[ref_label]]

# Align all to baseline using greedy matching
align_labels <- function(ref_vec, test_vec) {
  valid <- !is.na(ref_vec) & !is.na(test_vec)
  ref_v <- ref_vec[valid]
  test_v <- test_vec[valid]
  n_clust <- max(c(ref_v, test_v))
  if (n_clust > 6) return(1:n_clust)
  cm <- matrix(0, n_clust, n_clust)
  for (i in seq_along(ref_v)) cm[test_v[i], ref_v[i]] <- cm[test_v[i], ref_v[i]] + 1
  used_ref <- rep(FALSE, n_clust)
  mapping <- rep(NA, n_clust)
  for (iter in 1:n_clust) {
    best_val <- -1; best_t <- 0; best_r <- 0
    for (t in 1:n_clust) {
      if (!is.na(mapping[t])) next
      for (r in 1:n_clust) {
        if (used_ref[r]) next
        if (cm[t, r] > best_val) { best_val <- cm[t, r]; best_t <- t; best_r <- r }
      }
    }
    mapping[best_t] <- best_r; used_ref[best_r] <- TRUE
  }
  mapping
}

# Get sample names from FAM files
get_fam <- function(file) {
  read.table(paste0(file, ".fam"), header=FALSE)[,2]
}

# Build a unified data frame of all assignments (only samples in baseline)
baseline_samples <- get_fam(baseline_file)
n_baseline <- length(baseline_samples)

assign_df <- data.frame(Sample = baseline_samples, stringsAsFactors=FALSE)
assign_df$Baseline <- ref

for (i in 1:nrow(scenarios)) {
  lab <- scenarios$Label[i]
  if (lab == ref_label) next
  f <- scenarios$File[i]
  fam <- get_fam(f)
  test <- assignments[[lab]]
  names(test) <- fam
  
  # Find common samples with baseline
  common <- intersect(baseline_samples, fam)
  test_aligned <- rep(NA, n_baseline)
  test_aligned[match(common, baseline_samples)] <- test[common]
  
  # Align labels using common samples
  map <- align_labels(ref[match(common, baseline_samples)], test[common])
  test_aligned_valid <- !is.na(test_aligned)
  test_aligned[test_aligned_valid] <- map[test_aligned[test_aligned_valid]]
  
  assign_df[[lab]] <- test_aligned
}

# Cluster names (rice subpopulations)
pop_names <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")
# We'll try to infer population names from Q values later

# ============================================================
# 2. CONFUSION MATRIX (Analysis 1)
# ============================================================
cat("Computing confusion matrices...\n")

# For each scenario, create a confusion matrix vs baseline
confusion_list <- list()
change_summary <- data.frame(Scenario=character(), OverallAgreement=numeric(),
                             stringsAsFactors=FALSE)

for (lab in setdiff(scenarios$Label, ref_label)) {
  valid <- !is.na(assign_df$Baseline) & !is.na(assign_df[[lab]])
  if (sum(valid) == 0) next
  ct <- table(Baseline = assign_df$Baseline[valid], Scenario = assign_df[[lab]][valid])
  confusion_list[[lab]] <- ct
  change_summary <- rbind(change_summary, 
    data.frame(Scenario = lab, OverallAgreement = sum(diag(ct)) / sum(ct), stringsAsFactors=FALSE))
}

# Build per-accession change tracking
cat("Building per-accession change tracking...\n")
# For each baseline cluster, show how many accessions go to each new cluster per scenario
cluster_names <- paste0("C", 1:5)

change_tables <- list()
for (c in 1:5) {
  idx <- which(assign_df$Baseline == c)
  mat <- matrix(0, nrow=length(setdiff(scenarios$Label, ref_label)), ncol=5)
  rownames(mat) <- setdiff(scenarios$Label, ref_label)
  colnames(mat) <- cluster_names
  
  for (i in seq_along(setdiff(scenarios$Label, ref_label))) {
    lab <- setdiff(scenarios$Label, ref_label)[i]
    vals <- assign_df[[lab]][idx]
    for (j in 1:5) {
      mat[i, j] <- sum(vals == j, na.rm=TRUE)
    }
  }
  # Convert to proportions
  n_c <- length(idx)
  mat_prop <- mat / n_c * 100
  change_tables[[c]] <- list(count=mat, prop=mat_prop, n=n_c)
}

# ============================================================
# 3. PAIRWISE ARI AND NMI (Analysis 2)
# ============================================================
cat("Computing pairwise ARI and NMI...\n")

# For pairwise comparisons, use only samples common to both
pairwise_ari <- matrix(NA, nrow(scenarios), nrow(scenarios))
pairwise_nmi <- matrix(NA, nrow(scenarios), nrow(scenarios))
rownames(pairwise_ari) <- scenarios$Label
colnames(pairwise_ari) <- scenarios$Label
rownames(pairwise_nmi) <- scenarios$Label
colnames(pairwise_nmi) <- scenarios$Label

# Compute sample intersection for each scenario
sample_lists <- list()
for (i in 1:nrow(scenarios)) {
  f <- scenarios$File[i]
  sample_lists[[scenarios$Label[i]]] <- get_fam(f)
}

for (i in 1:nrow(scenarios)) {
  for (j in i:nrow(scenarios)) {
    lab_i <- scenarios$Label[i]
    lab_j <- scenarios$Label[j]
    
    common <- intersect(sample_lists[[lab_i]], sample_lists[[lab_j]])
    if (length(common) < 10) next
    
    # Get assignments for common samples
    fam_i <- get_fam(scenarios$File[i])
    fam_j <- get_fam(scenarios$File[j])
    
    ai <- assignments[[lab_i]][match(common, fam_i)]
    aj <- assignments[[lab_j]][match(common, fam_j)]
    
    # Align: map labels of j to i
    map <- align_labels(ai, aj)
    aj_aligned <- map[aj]
    
    # ARI
    tbl <- table(ai, aj_aligned)
    ari <- tryCatch({ 
      n <- sum(tbl)
      sum_i <- rowSums(tbl)
      sum_j <- colSums(tbl)
      n_choose_2 <- function(x) x*(x-1)/2
      index <- sum(sapply(tbl, n_choose_2))
      expected <- sum(sapply(sum_i, n_choose_2)) * sum(sapply(sum_j, n_choose_2)) / n_choose_2(n)
      max_index <- (sum(sapply(sum_i, n_choose_2)) + sum(sapply(sum_j, n_choose_2))) / 2
      (index - expected) / (max_index - expected)
    }, error=function(e) NA)
    
    # NMI
    nmi <- tryCatch({
      px <- rowSums(tbl) / sum(tbl)
      py <- colSums(tbl) / sum(tbl)
      pxy <- tbl / sum(tbl)
      mi <- sum(pxy * log(pxy / (outer(px, py, FUN=`*`)) + 1e-15), na.rm=TRUE)
      hx <- -sum(px * log(px + 1e-15), na.rm=TRUE)
      hy <- -sum(py * log(py + 1e-15), na.rm=TRUE)
      2 * mi / (hx + hy)
    }, error=function(e) NA)
    
    pairwise_ari[i,j] <- pairwise_ari[j,i] <- ari
    pairwise_nmi[i,j] <- pairwise_nmi[j,i] <- nmi
  }
}
diag(pairwise_ari) <- 1
diag(pairwise_nmi) <- 1

# ============================================================
# 4. PCA SIMILARITY (Analysis 3)
# ============================================================
cat("Computing PCA similarity...\n")

read_eigenvec <- function(file) {
  read.table(paste0(file, "_pca.eigenvec"), header=FALSE, stringsAsFactors=FALSE)
}

pca_cor <- matrix(NA, nrow(scenarios), nrow(scenarios))
rownames(pca_cor) <- scenarios$Label
colnames(pca_cor) <- scenarios$Label
pca_cor2 <- pca_cor

for (i in 1:nrow(scenarios)) {
  for (j in i:nrow(scenarios)) {
    lab_i <- scenarios$Label[i]
    lab_j <- scenarios$Label[j]
    f_i <- scenarios$File[i]
    f_j <- scenarios$File[j]
    
    pca_i <- read_eigenvec(f_i)
    pca_j <- read_eigenvec(f_j)
    
    # Match by sample ID
    id_i <- paste(pca_i$V1, pca_i$V2, sep="_")
    id_j <- paste(pca_j$V1, pca_j$V2, sep="_")
    common <- intersect(id_i, id_j)
    if (length(common) < 10) next
    
    idx_i <- match(common, id_i)
    idx_j <- match(common, id_j)
    
    pca_cor[i,j] <- pca_cor[j,i] <- cor(pca_i$V3[idx_i], pca_j$V3[idx_j], use="complete.obs")
    pca_cor2[i,j] <- pca_cor2[j,i] <- cor(pca_i$V4[idx_i], pca_j$V4[idx_j], use="complete.obs")
  }
}
diag(pca_cor) <- 1
diag(pca_cor2) <- 1

# ============================================================
# 5. GENETIC DIVERSITY: PIC, He, Ho, Shannon (Analysis 4)
# ============================================================
cat("Computing genetic diversity metrics...\n")

calc_pic <- function(maf) {
  # For biallelic markers: PIC = 1 - (p^2 + q^2) = 2*p*q = 2*MAF*(1-MAF)
  # Same as expected heterozygosity He
  2 * maf * (1 - maf)
}

calc_shannon <- function(maf) {
  p <- c(maf, 1-maf)
  -sum(p * log(p + 1e-15))
}

diversity_results <- data.frame(Scenario = scenarios$Label, Short = scenarios$Short,
                                Dim = scenarios$Dim, N_SNPs = NA,
                                Mean_MAF = NA, Mean_PIC = NA, Mean_He = NA,
                                Mean_Shannon = NA, Mean_Ho = NA,
                                stringsAsFactors=FALSE)

for (i in 1:nrow(scenarios)) {
  f <- scenarios$File[i]
  
  # Read frequency file
  frq <- tryCatch(read.table(paste0(f, "_freq.frq"), header=TRUE, stringsAsFactors=FALSE), 
                  error=function(e) NULL)
  if (!is.null(frq)) {
    maf <- ifelse(frq$MAF > 0.5, 1 - frq$MAF, frq$MAF)
    diversity_results$N_SNPs[i] <- nrow(frq)
    diversity_results$Mean_MAF[i] <- mean(maf, na.rm=TRUE)
    diversity_results$Mean_PIC[i] <- mean(calc_pic(maf), na.rm=TRUE)
    diversity_results$Mean_He[i] <- mean(calc_pic(maf), na.rm=TRUE)  # He = PIC for biallelic
    diversity_results$Mean_Shannon[i] <- mean(sapply(maf, calc_shannon), na.rm=TRUE)
  }
  
  # Read HWE for observed heterozygosity
  hwe <- tryCatch(read.table(paste0(f, "_hwe.hwe"), header=TRUE, stringsAsFactors=FALSE),
                  error=function(e) NULL)
  if (!is.null(hwe)) {
    # Ho = observed heterozygosity
    # PLINK HWE has O(HET) = observed heterozygotes
    n_total <- hwe$O.HET. + hwe$E.HET.  # approximate: O(HET) + O(HOM) = N
    # Actually, PLINK HWE reports: GENO, O(HOM A1), O(HET), E(HOM A1), E(HET), CHISQ, P
    # Let me check the column names
    # O(HET) = observed number of heterozygotes
    # We need total N per SNP to compute Ho = O(HET) / N
    # Use O(HOM A1) + O(HET) + O(HOM A2)... but PLINK only reports O(HOM A1) and O(HET)
    # Actually, the columns in PLINK .hwe are: CHR, SNP, TEST, A1, A2, GENO, O(HET), E(HET), P
    # where GENO is like "A1/A1, A1/A2, A2/A2" counts
    # We can parse GENO or use O(HET) / N
    # Simplest: total N = 2 * O(HET) / observed_het_freq... no.
    # Let me extract from GENO field instead
    geno_parts <- strsplit(hwe$GENO, "/")
    n1 <- as.numeric(sapply(geno_parts, `[`, 1))
    n_het <- as.numeric(sapply(geno_parts, `[`, 2))
    n2 <- as.numeric(sapply(geno_parts, `[`, 3))
    n_total <- n1 + n_het + n2
    ho <- n_het / n_total
    diversity_results$Mean_Ho[i] <- mean(ho, na.rm=TRUE)
  }
}

# ============================================================
# 6. FST COMPARISON (Analysis 5)
# ============================================================
cat("Analyzing FST...\n")

fst_data <- tryCatch(read.csv("fst_summary.csv", stringsAsFactors=FALSE), error=function(e) NULL)
if (!is.null(fst_data)) {
  # Map scenario labels
  scenario_map_fst <- c("L4_LD100"="L4", "L0_LDnone"="L0", "S3_Mind10Before"="S3",
                        "Gnone_GenoNone"="G0", "G4_Geno02"="G4", "M1_MAF01"="M1",
                        "G2_Geno10"="G2", "M0_MAFnone"="M0", "S4_Mind02"="S4",
                        "S0_NoMind"="S0", "L1_LD08"="L1", "S2_Mind05"="S2",
                        "M4_MAF10"="M4", "S1_Mind10After"="S1", "L2_LD05"="L2",
                        "M2_MAF03"="M2")
  
  fst_data$Short <- scenario_map_fst[fst_data$Scenario]
  fst_data$Baseline <- ifelse(fst_data$Short == "S3", "Baseline", "Scenario")
}

# ============================================================
# 7. AMOVA (Analysis 6)
# ============================================================
cat("Running AMOVA...\n")

# Need distance matrices. Check if they exist
amova_results <- data.frame(Scenario = character(), Short = character(),
                            Among_Pct = numeric(), Within_Pct = numeric(),
                            Phi_ST = numeric(), stringsAsFactors=FALSE)

run_amova <- function(file, label, short, pop_assignments) {
  dist_file <- paste0(file, "_dist.dist")
  if (!file.exists(dist_file)) return(NULL)
  
  # Read distance matrix
  dist_mat <- as.dist(as.matrix(read.table(dist_file, header=FALSE)))
  
  # Get population for each sample
  fam <- read.table(paste0(file, ".fam"), header=FALSE)
  common <- intersect(fam$V2, names(pop_assignments))
  if (length(common) < 10) return(NULL)
  
  idx <- match(common, fam$V2)
  pops <- factor(pop_assignments[match(common, names(pop_assignments))])
  
  # Subset distance matrix
  dist_sub <- as.dist(as.matrix(dist_mat)[idx, idx])
  
  # Run AMOVA using pegas
  amova_result <- pegas::amova(dist_sub ~ pops, nperm = 0)
  
  # Extract variance components  
  varcomp <- amova_result$varcomp
  among_pct <- varcomp[1] / sum(varcomp) * 100
  within_pct <- varcomp[2] / sum(varcomp) * 100
  phi_st <- varcomp[1] / sum(varcomp)
  
  data.frame(Scenario=label, Short=short, Among_Pct=among_pct, Within_Pct=within_pct, Phi_ST=phi_st,
             stringsAsFactors=FALSE)
}

# Use baseline cluster assignments as population labels
pop_labels <- ref
names(pop_labels) <- baseline_samples

# Check if distance files exist, if not, create them
for (i in 1:nrow(scenarios)) {
  f <- scenarios$File[i]
  lab <- scenarios$Label[i]
  short <- scenarios$Short[i]
  dist_file <- paste0(f, "_dist.dist")
  
  if (!file.exists(dist_file)) {
    cat("  Distance file missing for", lab, "- run PLINK first\n")
    next
  }
  
  res <- run_amova(f, lab, short, pop_labels)
  if (!is.null(res)) {
    amova_results <- rbind(amova_results, res)
  }
}

# ============================================================
# 8. GENERATE OUTPUT TABLES
# ============================================================
cat("Generating output tables...\n")

# Table: Confusion matrices summary
confusion_summary <- do.call(rbind, lapply(names(confusion_list), function(lab) {
  ct <- confusion_list[[lab]]
  data.frame(Scenario=lab, 
             rbind(ct),
             stringsAsFactors=FALSE)
}))

# Table: Change tracking per cluster
change_long <- do.call(rbind, lapply(1:5, function(c) {
  t <- change_tables[[c]]
  melt(t$prop, varnames=c("Scenario", "NewCluster"), value.name="Percent") %>%
    mutate(OriginalCluster = paste0("C", c), N = t$n)
}))
write.csv(change_long, "Table_S6_cluster_flow.csv", row.names=FALSE)

# Table: Pairwise ARI
write.csv(round(pairwise_ari, 4), "Table_S7_pairwise_ARI.csv")
write.csv(round(pairwise_nmi, 4), "Table_S8_pairwise_NMI.csv")

# Table: PCA correlations
pca_cor_df <- melt(round(pca_cor, 4), varnames=c("Scenario1", "Scenario2"), value.name="PC1_Corr")
pca_cor2_df <- melt(round(pca_cor2, 4), varnames=c("Scenario1", "Scenario2"), value.name="PC2_Corr")
pca_all <- merge(pca_cor_df, pca_cor2_df, by=c("Scenario1", "Scenario2"))
write.csv(pca_all, "Table_S9_PCA_correlations.csv", row.names=FALSE)

# Table: Genetic diversity
write.csv(diversity_results, "Table_S10_genetic_diversity.csv", row.names=FALSE)

# Table: FST comparison
if (!is.null(fst_data)) {
  write.csv(fst_data, "Table_S11_FST_comparison.csv", row.names=FALSE)
}

# Table: AMOVA
if (nrow(amova_results) > 0) {
  write.csv(amova_results, "Table_S12_AMOVA.csv", row.names=FALSE)
}

# ============================================================
# 9. FIGURES
# ============================================================
cat("Generating figures...\n")

# Figure 8: Confusion matrix heatmap - for each dimension
plot_confusion_heatmap <- function(dim_name, dim_label) {
  labs <- scenarios$Label[scenarios$Dim == dim_name & scenarios$Label != ref_label]
  if (length(labs) == 0) return(NULL)
  
  plots <- list()
  for (lab in labs) {
    ct <- confusion_list[[lab]]
    ct_long <- melt(ct, varnames=c("Baseline", "Scenario"))
    ct_long$Baseline <- factor(ct_long$Baseline, levels=1:5)
    ct_long$Scenario <- factor(ct_long$Scenario, levels=1:5)
    
    p <- ggplot(ct_long, aes(x=Scenario, y=Baseline, fill=value)) +
      geom_tile() +
      geom_text(aes(label=value), size=3) +
      scale_fill_gradient(low="white", high="steelblue") +
      labs(title=lab, x="Scenario Cluster", y="Baseline Cluster") +
      theme_minimal() +
      coord_fixed()
    plots[[lab]] <- p
  }
  grid.arrange(grobs=plots, ncol=min(length(plots), 3), top=dim_label)
}

pdf("Figure_8_confusion_matrix.pdf", width=12, height=8)
for (dim in unique(scenarios$Dim[scenarios$Dim != ""])) {
  labs <- scenarios$Label[scenarios$Dim == dim & scenarios$Label != ref_label]
  if (length(labs) == 0) next
  
  par(mfrow=c(2, min(3, ceiling(length(labs)/2))))
  for (lab in labs) {
    ct <- confusion_list[[lab]]
    # Convert to percentage
    ct_pct <- sweep(ct, 1, rowSums(ct), "/") * 100
    
    # Color heatmap
    image(1:ncol(ct), 1:nrow(ct), t(ct_pct[nrow(ct):1,]), 
          col=colorRampPalette(c("white", "steelblue", "navy"))(100),
          xlab="Scenario Cluster", ylab="Baseline Cluster",
          main=lab, axes=FALSE)
    axis(1, at=1:ncol(ct), labels=1:ncol(ct))
    axis(2, at=1:nrow(ct), labels=nrow(ct):1)
    
    # Add text
    for (i in 1:nrow(ct)) {
      for (j in 1:ncol(ct)) {
        text(j, nrow(ct)-i+1, round(ct_pct[i,j], 1), cex=0.8, col=ifelse(ct_pct[i,j] > 50, "white", "black"))
      }
    }
    mtext(paste0("Agreement: ", round(sum(diag(ct))/sum(ct)*100, 1), "%"), side=3, cex=0.8)
  }
}
dev.off()
cat("Figure 8 (confusion matrix) saved\n")

# Figure 9: ARI + NMI heatmaps
pdf("Figure_9_ARI_NMI.pdf", width=14, height=6)
par(mfrow=c(1,2))

# ARI heatmap
ari_melt <- melt(pairwise_ari, varnames=c("S1", "S2"), value.name="ARI")
ari_melt$S1 <- factor(ari_melt$S1, levels=rev(scenarios$Label))
ari_melt$S2 <- factor(ari_melt$S2, levels=scenarios$Label)

p1 <- ggplot(ari_melt, aes(x=S2, y=S1, fill=ARI)) +
  geom_tile() +
  geom_text(aes(label=round(ARI, 2)), size=2.5) +
  scale_fill_gradient(low="white", high="darkgreen", limits=c(0,1), na.value="gray90") +
  labs(title="Adjusted Rand Index", x="", y="") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7),
        axis.text.y=element_text(size=7))

# NMI heatmap
nmi_melt <- melt(pairwise_nmi, varnames=c("S1", "S2"), value.name="NMI")
nmi_melt$S1 <- factor(nmi_melt$S1, levels=rev(scenarios$Label))
nmi_melt$S2 <- factor(nmi_melt$S2, levels=scenarios$Label)

p2 <- ggplot(nmi_melt, aes(x=S2, y=S1, fill=NMI)) +
  geom_tile() +
  geom_text(aes(label=round(NMI, 2)), size=2.5) +
  scale_fill_gradient(low="white", high="darkorange", limits=c(0,1), na.value="gray90") +
  labs(title="Normalized Mutual Information", x="", y="") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7),
        axis.text.y=element_text(size=7))

grid.arrange(p1, p2, ncol=2)
dev.off()
cat("Figure 9 (ARI/NMI) saved\n")

# Figure 10: PCA correlation heatmap
pdf("Figure_10_PCA_correlations.pdf", width=10, height=8)
par(mfrow=c(1,2))

# PC1
pca1_melt <- melt(pca_cor, varnames=c("S1", "S2"), value.name="Cor")
image(1:nrow(pca_cor), 1:ncol(pca_cor), pca_cor[nrow(pca_cor):1,],
      col=colorRampPalette(c("white", "orange", "darkred"))(100),
      xlab="", ylab="", main="PC1 Correlation", axes=FALSE)
axis(1, at=1:ncol(pca_cor), labels=colnames(pca_cor), las=2, cex.axis=0.5)
axis(2, at=1:nrow(pca_cor), labels=rownames(pca_cor)[nrow(pca_cor):1], las=2, cex.axis=0.5)

# PC2
image(1:nrow(pca_cor2), 1:ncol(pca_cor2), pca_cor2[nrow(pca_cor2):1,],
      col=colorRampPalette(c("white", "orange", "darkred"))(100),
      xlab="", ylab="", main="PC2 Correlation", axes=FALSE)
axis(1, at=1:ncol(pca_cor2), labels=colnames(pca_cor2), las=2, cex.axis=0.5)
axis(2, at=1:nrow(pca_cor2), labels=rownames(pca_cor2)[nrow(pca_cor2):1], las=2, cex.axis=0.5)

dev.off()
cat("Figure 10 (PCA correlations) saved\n")

# Figure 11: Genetic diversity comparison
pdf("Figure_11_genetic_diversity.pdf", width=12, height=8)
par(mfrow=c(2,3))
baseline_val <- diversity_results[diversity_results$Short == "S3", ]
for (metric in c("Mean_PIC", "Mean_He", "Mean_Ho", "Mean_Shannon", "Mean_MAF", "N_SNPs")) {
  vals <- diversity_results[[metric]][diversity_results$Short == "S3"]
  bl <- ifelse(length(vals) > 0, vals[1], NA)
  
  for (dim in c("Sample", "MAF", "GENO", "LD")) {
    sub <- diversity_results[diversity_results$Dim == dim, ]
    if (nrow(sub) <= 1) next
    
    bp <- barplot(sub[[metric]], names.arg=sub$Short, las=2, cex.names=0.8,
                  main=paste(metric, "-", dim), ylab=metric, cex.main=0.9,
                  col=ifelse(sub$Short == "S3", "red", "steelblue"))
    if (!is.na(bl)) abline(h=bl, lty=2, col="red", lwd=2)
  }
}
dev.off()
cat("Figure 11 (genetic diversity) saved\n")

# Figure 12: FST comparison
if (!is.null(fst_data)) {
  pdf("Figure_12_FST_comparison.pdf", width=10, height=6)
  
  par(mfrow=c(1,2))
  
  # Barplot by dimension
  dim_colors <- c(Sample="skyblue", MAF="salmon", GENO="lightgreen", LD="plum")
  fst_data$Dim <- ifelse(grepl("^S", fst_data$Short), "Sample",
                         ifelse(grepl("^M", fst_data$Short), "MAF",
                                ifelse(grepl("^G", fst_data$Short), "GENO", "LD")))
  
  # Sort by dimension then FST
  fst_data <- fst_data[order(fst_data$Dim, fst_data$Fst_mean), ]
  
  barplot(fst_data$Fst_mean, names.arg=fst_data$Short, las=2, cex.names=0.7,
          col=dim_colors[fst_data$Dim], main="Mean FST by Scenario",
          ylab="Mean FST", ylim=c(0, max(fst_data$Fst_mean)*1.2))
  abline(h=fst_data$Fst_mean[fst_data$Short == "S3"], lty=2, col="red", lwd=2)
  legend("topright", fill=dim_colors, legend=names(dim_colors), cex=0.7)
  
  # Line plot faceted by dimension
  plot(fst_data$Fst_mean, type="n", xlab="Scenario", ylab="Mean FST",
       main="FST Relative to Baseline", xaxt="n")
  for (dim in unique(fst_data$Dim)) {
    sub <- fst_data[fst_data$Dim == dim, ]
    lines(seq_along(sub$Short), sub$Fst_mean, type="o", pch=19, 
          col=dim_colors[dim], lwd=2)
  }
  axis(1, at=1:nrow(fst_data), labels=fst_data$Short, las=2, cex.axis=0.7)
  abline(h=fst_data$Fst_mean[fst_data$Short == "S3"], lty=2, col="gray", lwd=2)
  legend("topleft", fill=dim_colors, legend=names(dim_colors), cex=0.7)
  
  dev.off()
  cat("Figure 12 (FST comparison) saved\n")
}

# Figure 13: AMOVA
if (nrow(amova_results) > 0) {
  pdf("Figure_13_AMOVA.pdf", width=8, height=6)
  
  dim_colors <- c(Sample="skyblue", MAF="salmon", GENO="lightgreen", LD="plum")
  amova_results$Dim <- ifelse(grepl("^S", amova_results$Short), "Sample",
                              ifelse(grepl("^M", amova_results$Short), "MAF",
                                     ifelse(grepl("^G", amova_results$Short), "GENO", "LD")))
  
  barplot(amova_results$Among_Pct, names.arg=amova_results$Short, las=2, cex.names=0.7,
          col=dim_colors[amova_results$Dim], 
          main="Among-population Variance (%)", ylab="% Variance",
          ylim=c(0, 100))
  abline(h=amova_results$Among_Pct[amova_results$Short == "S3"], lty=2, col="red", lwd=2)
  legend("topright", fill=dim_colors, legend=names(dim_colors), cex=0.7)
  
  dev.off()
  cat("Figure 13 (AMOVA) saved\n")
}

# ============================================================
# 10. CLUSTER FLOW FIGURE (Supplement to Confusion Matrix)
# ============================================================
# Show how each baseline cluster redistributes across scenarios
pdf("Figure_S5_cluster_flow.pdf", width=12, height=8)
par(mfrow=c(2,3))
for (c in 1:5) {
  t <- change_tables[[c]]
  d <- as.data.frame(t$prop)
  d$Scenario <- rownames(d)
  d_long <- melt(d, id.vars="Scenario", variable.name="NewCluster", value.name="Percent")
  
  barplot(t(t$prop), beside=TRUE, col=rainbow(5), 
          main=paste("Baseline Cluster", c, "(n=", t$n, ")"),
          xlab="Scenario", ylab="% reassigned", las=2, cex.names=0.7)
  legend("topright", fill=rainbow(5), legend=paste("C", 1:5), cex=0.6)
}
dev.off()
cat("Figure S5 (cluster flow) saved\n")

cat("\n===== ALL ANALYSES COMPLETE =====\n")
cat("Tables: S6 (cluster flow), S7 (ARI), S8 (NMI), S9 (PCA cor), S10 (diversity), S11 (FST), S12 (AMOVA)\n")
cat("Figures: 8 (confusion), 9 (ARI/NMI), 10 (PCA cor), 11 (diversity), 12 (FST), 13 (AMOVA), S5 (flow)\n")
