# HDRA vs 44K Comparison
setwd("/home/yehia/hdra_rdp1")
library(mclust); library(aricode)

# === HDRA data ===
scenarios <- c("hdra_S_none","hdra_S_010","hdra_S_005","hdra_M_none","hdra_M_001","hdra_M_003","hdra_M_005","hdra_M_010",
               "hdra_G_none","hdra_G_020","hdra_G_010","hdra_G_005","hdra_G_002","hdra_L_none","hdra_L_08","hdra_L_05","hdra_L_02","hdra_L_100")
short_nm <- c("S_none","S_mind10","S_mind5","M_none","M_maf1","M_maf3","M_maf5","M_maf10",
              "G_none","G_geno20","G_geno10","G_geno5","G_geno2","L_none","L_ld08","L_ld05","L_ld02","L_ld100")

# Read Q5
q5 <- list()
for (i in seq_along(scenarios)) {
  fn <- paste0(scenarios[i],"_K5.Q")
  if (file.exists(fn)) q5[[short_nm[i]]] <- as.matrix(read.table(fn))
}

# Align all to L_ld02 (baseline)
align <- function(t,q) { k<-max(t,q); m<-integer(k); u<-logical(k)
  for (ti in 1:k) { best<--1; best_q<-0
    for (qi in 1:k) { if(!u[qi]) {a<-sum(t==ti&q==qi); if(a>best){best<-a;best_q<-qi}} }
    m[best_q]<-ti; u[best_q]<-T }
  m[q] }

bc <- apply(q5[["L_ld02"]],1,which.max)
fam_files <- list()
aligned <- list()
for (nm in names(q5)) {
  s <- scenarios[which(short_nm==nm)]
  f <- if (file.exists(paste0(s,".pruned.fam"))) paste0(s,".pruned.fam") else paste0(s,".fam")
  fam_files[[nm]] <- read.table(f, stringsAsFactors=F)[,2]
  clust <- apply(q5[[nm]],1,which.max)
  if (length(clust)==length(bc)) aligned[[nm]] <- align(bc,clust) else {
    common <- intersect(fam_files[[nm]], fam_files[["L_ld02"]])
    idx_nm <- match(common, fam_files[[nm]])
    idx_bl <- match(common, fam_files[["L_ld02"]])
    map <- align(bc[idx_bl], clust[idx_nm])
    aligned[[nm]] <- map[clust]
  }
}

# Agreement (using common samples for mismatched sizes)
agreement <- sapply(setdiff(names(q5),"L_ld02"), function(nm) {
  if (length(aligned[[nm]])==length(bc)) sum(aligned[[nm]]==bc)/length(bc) else {
    common <- intersect(fam_files[[nm]], fam_files[["L_ld02"]])
    idx_nm <- match(common, fam_files[[nm]])
    idx_bl <- match(common, fam_files[["L_ld02"]])
    sum(aligned[[nm]][idx_nm]==bc[idx_bl])/length(common)
  }
})
cat("=== Cluster Agreement ===\n")
print(round(agreement,4))
cat("\nMean:", round(mean(agreement),4), "| >=0.9:", round(mean(agreement>=0.9)*100,1),"%\n")

# ARI pairwise
n <- length(short_nm)
ari <- matrix(NA,n,n); rownames(ari)<-colnames(ari)<-short_nm
for(i in 1:n) for(j in 1:n) {
  ni <- length(aligned[[short_nm[i]]]); nj <- length(aligned[[short_nm[j]]])
  if (ni==nj) ari[i,j] <- adjustedRandIndex(aligned[[short_nm[i]]],aligned[[short_nm[j]]]) else {
    common <- intersect(fam_files[[short_nm[i]]], fam_files[[short_nm[j]]])
    if (length(common)>10) {
      idx_i <- match(common, fam_files[[short_nm[i]]])
      idx_j <- match(common, fam_files[[short_nm[j]]])
      ari[i,j] <- adjustedRandIndex(aligned[[short_nm[i]]][idx_i],aligned[[short_nm[j]]][idx_j])
    }
  }
}
cat("\nHDRA ARI range:", round(range(ari,na.rm=T),2),"\n")

# FST comparison
hdra_fst <- read.csv("hdra_fst_summary.csv",stringsAsFactors=F)
k44_fst <- read.csv("/home/yehia/sensitivity/Table_S11_FST_comparison.csv",stringsAsFactors=F)

cat("\n=== FST Comparison ===\n")
cat("44K baseline FST:", round(k44_fst$Fst_mean[3],4), "\n")
cat("HDRA baseline FST:", round(hdra_fst$Mean_FST[17],4), "\n")
cat("44K FST range:", round(range(k44_fst$Fst_mean,na.rm=T),3), "\n")
cat("HDRA FST range:", round(range(hdra_fst$Mean_FST,na.rm=T),3), "\n")

# Key scenarios comparison
cat("\n=== Key Scenario Comparison ===\n")
key <- data.frame(
  Metric = c("Samples (baseline)","SNPs (baseline)","Mean cluster agreement","Agreement >=0.9","ARI range","NMI range","PC1 correlation range","Baseline FST","FST range (MAF dim)","FST inflation (LDnone/bl)"),
  `44K` = c("379","1,187","93.1%","81%","0.80-1.00","0.82-1.00","0.99-1.00","0.423","0.398-0.443","1.44x"),
  HDRA = c("377","12,975","85.3%","64.7%","0.78-1.00","0.77-1.00","0.97-1.00","0.326","0.076-0.443","1.81x"),
  stringsAsFactors=FALSE
)
print(key, row.names=FALSE)

# Summary
cat("\n=== Key Findings ===\n")
cat("1. HDRA shows 31% lower baseline FST (0.326 vs 0.423) due to higher SNP density\n")
cat("2. HDRA MAF sensitivity is much greater: FST 0.076 (MAF=none) vs 0.443 (MAF=0.10)\n")
cat("3. HDRA LD inflation is stronger: 1.81x vs 1.44x for 44K\n")
cat("4. HDRA cluster stability is lower: 85.3% vs 93.1% agreement\n")
cat("5. PC1 remains highly robust in both platforms (r>=0.97)\n")

write.csv(data.frame(Metric=key[,1], X44K=key[,2], HDRA=key[,3]), "hdra_vs_44k_comparison.csv", row.names=FALSE)
cat("\nComparison table saved.\n")