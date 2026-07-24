# HDRA vs 44K Comparison
library(mclust); library(aricode)

setwd("/home/yehia/hdra_rdp1")

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

# Read fam files
read_fam <- function(nm) {
  s <- scenarios[which(short_nm==nm)]
  f <- paste0(s,".pruned.fam")
  if (!file.exists(f)) f <- paste0(s,".fam")
  if (!file.exists(f)) return(NULL)
  read.table(f, stringsAsFactors=F)[,2]
}
fams <- lapply(short_nm, read_fam)
names(fams) <- short_nm

# Align clusters
align_fn <- function(t, q) {
  if (length(t)==0 || length(q)==0) return(q)
  k <- max(t,q)
  m <- integer(k); u <- logical(k)
  for (ti in 1:k) {
    best <- -1; best_q <- 0
    for (qi in 1:k) {
      if (!u[qi]) { a <- sum(t==ti & q==qi); if (a>best) { best<-a; best_q<-qi } }
    }
    m[best_q] <- ti; u[best_q] <- TRUE
  }
  m[q]
}

bc <- apply(q5[["L_ld02"]], 1, which.max)
aligned <- list()
for (nm in names(q5)) {
  clust <- apply(q5[[nm]], 1, which.max)
  if (length(clust)==length(bc)) {
    aligned[[nm]] <- align_fn(bc, clust)
  } else {
    common <- intersect(fams[[nm]], fams[["L_ld02"]])
    idx_nm <- match(common, fams[[nm]])
    idx_bl <- match(common, fams[["L_ld02"]])
    if (length(common) > 5) {
      map <- align_fn(bc[idx_bl], clust[idx_nm])
      aligned[[nm]] <- map[clust]
    } else {
      aligned[[nm]] <- clust
    }
  }
}

# Agreement
cat("=== Cluster Agreement ===\n")
agree <- sapply(setdiff(names(aligned),"L_ld02"), function(nm) {
  if (length(aligned[[nm]])==length(bc)) {
    sum(aligned[[nm]]==bc)/length(bc)
  } else {
    common <- intersect(fams[[nm]], fams[["L_ld02"]])
    idx_nm <- match(common, fams[[nm]])
    idx_bl <- match(common, fams[["L_ld02"]])
    sum(aligned[[nm]][idx_nm]==bc[idx_bl])/length(common)
  }
})
print(round(agree,4))
cat("Mean:", round(mean(agree),4), "| >=0.9:", round(mean(agree>=0.9)*100,1), "%\n")

# ARI
cat("\n=== Pairwise ARI ===\n")
n <- length(short_nm)
ari <- matrix(NA,n,n)
rownames(ari) <- colnames(ari) <- short_nm
for (i in 1:n) for (j in 1:n) {
  ni <- length(aligned[[short_nm[i]]])
  nj <- length(aligned[[short_nm[j]]])
  if (ni==nj && ni>0) {
    ari[i,j] <- adjustedRandIndex(aligned[[short_nm[i]]], aligned[[short_nm[j]]])
  } else if (ni>0 && nj>0) {
    common <- intersect(fams[[short_nm[i]]], fams[[short_nm[j]]])
    if (length(common)>5) {
      id_i <- match(common, fams[[short_nm[i]]])
      id_j <- match(common, fams[[short_nm[j]]])
      ari[i,j] <- adjustedRandIndex(aligned[[short_nm[i]]][id_i], aligned[[short_nm[j]]][id_j])
    }
  }
}
write.csv(ari, "hdra_ari.csv")
cat("ARI range:", round(range(ari,na.rm=T),2),"\n")

# FST
hdra_fst <- read.csv("hdra_fst_summary.csv", stringsAsFactors=F)
k44_fst <- read.csv("/home/yehia/sensitivity/Table_S11_FST_comparison.csv", stringsAsFactors=F)

cat("\n=== FST Comparison ===\n")
cat("44K baseline FST:", round(k44_fst$Fst_mean[3],4), "\n")
cat("HDRA baseline FST:", round(hdra_fst$Mean_FST[17],4), "\n")
cat("44K FST range:", round(range(k44_fst$Fst_mean,na.rm=T),3), "\n")
cat("HDRA FST range:", round(range(hdra_fst$Mean_FST,na.rm=T),3), "\n")

# LD inflation
bl44 <- k44_fst$Fst_mean[3]
ldnone44 <- k44_fst$Fst_mean[2]
bl_hdra <- hdra_fst$Mean_FST[17]
ldnone_hdra <- hdra_fst$Mean_FST[14]
cat("44K LD inflation (L_none/baseline):", round(ldnone44/bl44,2), "x\n")
cat("HDRA LD inflation:", round(ldnone_hdra/bl_hdra,2), "x\n")

# MAF range comparison
cat("\n44K MAF FST range:", round(range(k44_fst$Fst_mean[4:8],na.rm=T),3), "\n")
cat("HDRA MAF FST range:", round(range(hdra_fst$Mean_FST[4:8],na.rm=T),3), "\n")

# Summary table
cat("\n\n=== Key Metrics Comparison ===\n")
tbl <- data.frame(
  Metric = c("Baseline samples","Baseline SNPs","Mean cluster agreement","Proportion >=0.9","ARI range","Baseline FST","FST range (MAF dim)","LD inflation factor","PC1 correlation range"),
  X44K = c("379","1,187","93.1%","81%","0.80-1.00","0.423","0.398-0.443","1.44x","0.99-1.00"),
  HDRA = c("377","12,975","85.3%","64.7%","0.78-1.00","0.326","0.076-0.443","1.81x","0.97-1.00"),
  stringsAsFactors=F
)
print(tbl, row.names=F)

write.csv(tbl, "hdra_vs_44k_comparison.csv", row.names=F)
cat("\nTable saved.\n")
