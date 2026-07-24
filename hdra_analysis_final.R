# HDRA RDP1 Sensitivity Analysis - Fixed
setwd("/home/yehia/hdra_rdp1")
library(mclust)
library(aricode)

scenarios <- c(
  "hdra_S_none","hdra_S_010","hdra_S_005",
  "hdra_M_none","hdra_M_001","hdra_M_003","hdra_M_005","hdra_M_010",
  "hdra_G_none","hdra_G_020","hdra_G_010","hdra_G_005","hdra_G_002",
  "hdra_L_none","hdra_L_08","hdra_L_05","hdra_L_02","hdra_L_100"
)
short_names <- c(
  "S_none","S_mind10","S_mind5",
  "M_none","M_maf1","M_maf3","M_maf5","M_maf10",
  "G_none","G_geno20","G_geno10","G_geno5","G_geno2",
  "L_none","L_ld08","L_ld05","L_ld02","L_ld100"
)
dimensions <- c(rep("Sample\nfiltering",3), rep("MAF\nthreshold",5),
                rep("GENO\nthreshold",5), rep("LD\npruning",5))

# Read K=5 Q-files
q5 <- list()
for (i in seq_along(scenarios)) {
  fn <- paste0(scenarios[i], "_K5.Q")
  if (file.exists(fn)) q5[[short_names[i]]] <- as.matrix(read.table(fn))
}
baseline <- "L_ld02"

# Greedy cluster alignment
align <- function(t, q) {
  k <- max(t, q)
  m <- integer(k)
  u <- logical(k)
  for (ti in 1:k) {
    best <- -1; best_q <- 0
    for (qi in 1:k) { if (!u[qi]) { a <- sum(t==ti & q==qi); if (a>best) { best<-a; best_q<-qi } } }
    m[best_q] <- ti; u[best_q] <- TRUE
  }
  m[q]
}

aligned <- list()
bc <- apply(q5[[baseline]], 1, which.max)
for (nm in names(q5)) {
  q <- q5[[nm]]; clust <- apply(q, 1, which.max)
  aligned[[nm]] <- if (length(clust)==length(bc)) align(bc, clust) else clust
}

# Agreement
agreement <- sapply(setdiff(names(aligned), baseline), function(nm) {
  sum(aligned[[nm]] == aligned[[baseline]]) / length(aligned[[baseline]])
})
cat("Mean agreement:", round(mean(agreement),4), "\n")
cat("Agreement >=0.9:", round(mean(agreement>=0.9)*100,1), "%\n\n")

# ARI / NMI pairwise
n <- length(short_names)
ari <- matrix(NA,n,n); nmi <- ari
rownames(ari)<-colnames(ari)<-short_names; rownames(nmi)<-colnames(nmi)<-short_names
for (i in 1:n) for (j in 1:n) {
  if (length(aligned[[short_names[i]]])==length(aligned[[short_names[j]]])) {
    ari[i,j] <- adjustedRandIndex(aligned[[short_names[i]]], aligned[[short_names[j]]])
    nmi[i,j] <- NMI(aligned[[short_names[i]]], aligned[[short_names[j]]])
  }
}
write.csv(ari,"hdra_ari.csv"); write.csv(nmi,"hdra_nmi.csv")
cat("ARI range:", round(range(ari,na.rm=T),3),"\n")
cat("NMI range:", round(range(nmi,na.rm=T),3),"\n\n")

# PCA correlations
pca_cor <- matrix(NA,length(scenarios),length(scenarios))
rownames(pca_cor)<-colnames(pca_cor)<-short_names
for (i in seq_along(scenarios)) {
  f1 <- paste0(scenarios[i],"_pca.eigenvec")
  if (!file.exists(f1)) next
  p1 <- try(read.table(f1,header=F)[,3]); if(inherits(p1,"try-error")) next
  for (j in seq_along(scenarios)) {
    f2<-paste0(scenarios[j],"_pca.eigenvec")
    if(!file.exists(f2)) next
    p2<-try(read.table(f2,header=F)[,3]); if(inherits(p2,"try-error")) next
    if(length(p1)==length(p2)) pca_cor[i,j]<-abs(cor(p1,p2))
  }
}
write.csv(pca_cor,"hdra_pca_cor.csv")
cat("PC1 correlation range:", round(range(pca_cor,na.rm=T),3),"\n\n")

# Diversity from freq files (MAF -> He, PIC)
div_list <- list()
for (i in seq_along(scenarios)) {
  fn <- paste0(scenarios[i],"_freq.frq")
  if (!file.exists(fn)) next
  frq <- try(read.table(fn, header=TRUE, stringsAsFactors=FALSE))
  if (inherits(frq,"try-error") || nrow(frq)<10) next
  maf <- frq$MAF
  he <- mean(2*maf*(1-maf), na.rm=T)
  ho <- NA  # not in freq file
  pic <- mean(1 - (maf^2 + (1-maf)^2), na.rm=T)
  div_list[[short_names[i]]] <- c(He=round(he,4), PIC=round(pic,4))
}
div_df <- do.call(rbind, div_list)
write.csv(div_df,"hdra_diversity.csv")
cat("Diversity (He):\n"); print(round(div_df[,1],4))

# FST summary
if (file.exists("hdra_fst_summary.csv")) {
  fst <- read.csv("hdra_fst_summary.csv", stringsAsFactors=FALSE)
  cat("\nFST range:", round(range(fst$Mean_FST,na.rm=T),3),"\n")
  cat("Baseline FST (L_ld02):", round(fst$Mean_FST[which(short_names=="L_ld02")],4),"\n\n")
}

# Summary table
stats <- read.csv("s.csv", stringsAsFactors=FALSE)
available <- short_names %in% names(agreement)
valid_scenarios <- scenarios[available]
valid_short <- short_names[available]
valid_dim <- dimensions[available]

summary_df <- data.frame(
  Scenario = valid_short,
  Dimension = valid_dim,
  Samples = stats$samples[match(valid_scenarios, stats$scenario)],
  SNPs = stats$snps[match(valid_scenarios, stats$scenario)],
  Agreement = round(agreement[valid_short],4),
  He = div_df[valid_short, "He"],
  PIC = div_df[valid_short, "PIC"],
  row.names = NULL
)
write.csv(summary_df, "hdra_summary.csv", row.names=FALSE)
cat("Summary table saved.\n")
print(summary_df[,-c(6,7)])

# Compare with 44K
cat("\n\n=== Comparison with 44K ===\n")
if (file.exists("/home/yehia/sensitivity/summary.csv")) {
  k44 <- read.csv("/home/yehia/sensitivity/summary.csv")
  hidx <- which(valid_short=="L_ld02")
  cat("44K baseline SNPs: 1,187 | HDRA baseline SNPs:", summary_df$SNPs[hidx], "\n")
  cat("44K mean agreement: 93.1% | HDRA:", round(mean(agreement,na.rm=T)*100,1),"%\n")
  cat("44K ARI range: 0.80-1.00 | HDRA ARI range:", round(range(ari,na.rm=T),2),"\n")
}

cat("\n=== Done ===\n")
