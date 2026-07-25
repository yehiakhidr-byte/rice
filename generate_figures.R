# Generate all figures for cross-platform manuscript
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gridExtra)

dir.create("figures", showWarnings=FALSE)
cols <- c("44K"="#2166AC", "HDRA"="#B2182B")

# ---- Figure 2: SNP/Sample retention ----
ld44 <- data.frame(Pruning=c("None","r2<0.8","r2<0.5","50/5/0.2","100/10/0.2"),
                   SNPs=c(26474,9609,4278,1187,822))
ld44$Pruning <- factor(ld44$Pruning, levels=unique(ld44$Pruning))
p2a <- ggplot(ld44, aes(x=Pruning, y=SNPs)) +
  geom_bar(stat="identity", fill="#2166AC", alpha=0.8) +
  geom_text(aes(label=SNPs), vjust=-0.3, size=3) +
  labs(title="44K: SNP retention by LD pruning", y="SNPs retained") +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1))

hdra_fst <- read.csv("hdra_fst_summary.csv")
maf <- hdra_fst[4:8,]
maf$Lab <- c("None","0.01","0.03","0.05","0.10")
maf$Lab <- factor(maf$Lab, levels=unique(maf$Lab))
p2b <- ggplot(maf, aes(x=Lab, y=SNPs)) +
  geom_bar(stat="identity", fill="#B2182B", alpha=0.8) +
  geom_text(aes(label=SNPs), vjust=-0.3, size=3) +
  labs(title="HDRA: SNP retention by MAF", y="SNPs retained") +
  theme_bw()

pdf("figures/Figure2.pdf", width=10, height=5)
grid.arrange(p2a, p2b, ncol=2)
dev.off()
cat("Figure 2 saved\n")

# ---- Figure 3: Cross-platform comparison ----
# 3a: FST
k44 <- read.csv("/home/yehia/sensitivity/Table_S11_FST_comparison.csv")
idx_44 <- c(3,7,8,9,2,4,1,5)
idx_hd <- c(17,4,8,9,14,16,1,3)
labs <- c("Baseline","MAF none","MAF 0.10","GENO none","LD none","LD r2 0.5","Sample none","mind 0.05")
fst_df <- rbind(
  data.frame(Scenario=labs, FST=k44$Fst_mean[idx_44], Platform="44K"),
  data.frame(Scenario=labs, FST=hdra_fst$Mean_FST[idx_hd], Platform="HDRA")
)
fst_df$Scenario <- factor(fst_df$Scenario, levels=rev(unique(fst_df$Scenario)))

p3a <- ggplot(fst_df, aes(x=FST, y=Scenario, color=Platform)) +
  geom_point(size=3, position=position_dodge(width=0.5)) +
  scale_color_manual(values=cols) +
  labs(title="a  FST across scenarios", x="Mean FST") +
  theme_bw() + theme(legend.position="bottom")

# 3b: Cluster agreement HDRA
sc <- c("S_none","S_mind5","M_none","M_maf1","M_maf3","M_maf10","G_none","G_geno20","G_geno10","G_geno2","L_none","L_ld08","L_ld05","L_ld100")
dim <- c("Sample","Sample","MAF","MAF","MAF","MAF","GENO","GENO","GENO","GENO","LD","LD","LD","LD")
ag <- c(0.593,0.511,0.850,0.861,0.944,0.989,0.826,0.826,0.817,0.974,0.904,0.910,0.938,0.982)
agree_df <- data.frame(Scenario=sc, Agreement=ag, Dimension=dim)

p3b <- ggplot(agree_df, aes(x=Agreement, y=Scenario)) +
  geom_point(size=2.5, color="#B2182B") +
  facet_wrap(~Dimension, scales="free_y", ncol=1) +
  labs(title="b  HDRA cluster agreement", x="Proportion agreement") +
  xlim(0.4, 1.05) + theme_bw()

# 3c: Diversity
hdra_div <- read.csv("hdra_diversity.csv", row.names=1)
div_df <- data.frame(
  Lab=c("MAF none","MAF 0.01","MAF 0.03","MAF 0.05","MAF 0.10"),
  HDRA=hdra_div$He[4:8],
  X44K=c(0.30,0.30,0.30,0.30,0.27)
)
div_m <- melt(div_df, id.vars="Lab", variable.name="Platform", value.name="He")

p3c <- ggplot(div_m, aes(x=He, y=Lab, color=Platform)) +
  geom_point(size=3, position=position_dodge(width=0.5)) +
  scale_color_manual(values=cols) +
  labs(title="c  He by MAF threshold", x="Expected heterozygosity") +
  theme_bw() + theme(legend.position="bottom")

pdf("figures/Figure3.pdf", width=12, height=8)
grid.arrange(p3a, p3b, p3c, ncol=3, widths=c(1,1,1))
dev.off()
cat("Figure 3 saved\n")

# ---- Figure 4: PCA robustness ----
hdra_pca <- read.csv("hdra_pca_cor.csv", row.names=1)
bl_col <- grep("L_ld02", colnames(hdra_pca))[1]
pca_hdra <- data.frame(
  Scenario=rownames(hdra_pca),
  PC1=hdra_pca[,bl_col],
  Platform="HDRA"
)
pca_hdra <- pca_hdra[!is.na(pca_hdra$PC1),]
pca_hdra <- pca_hdra[pca_hdra$Scenario != "L_ld02",]

pca_44k <- data.frame(
  Scenario=c("Baseline","MAF none","MAF 0.10","GENO none","GENO 0.02","LD none","Sample none","mind 0.05"),
  PC1=c(1.00,0.995,0.995,0.995,0.995,0.994,0.999,0.82),
  Platform="44K"
)
pca_all <- rbind(pca_hdra, pca_44k)

p4 <- ggplot(pca_all, aes(x=PC1, y=Scenario, color=Platform)) +
  geom_point(size=3, position=position_dodge(width=0.5)) +
  scale_color_manual(values=cols) +
  geom_vline(xintercept=0.97, linetype="dashed", alpha=0.4) +
  labs(x="PC1 correlation with baseline", y="") +
  xlim(0.8, 1.02) + theme_bw() + theme(legend.position="bottom")

pdf("figures/Figure4.pdf", width=8, height=6)
print(p4)
dev.off()
cat("Figure 4 saved\n")

# ---- Figure 5: FST vs MAF ----
maf5 <- hdra_fst[4:8,]
maf5$MAFthresh <- c(0,0.01,0.03,0.05,0.10)
p5 <- ggplot(maf5, aes(x=MAFthresh, y=Mean_FST)) +
  geom_line(color="#B2182B", size=1) +
  geom_point(color="#B2182B", size=3) +
  geom_text(aes(label=round(Mean_FST,3)), vjust=-1, size=3.5) +
  scale_x_continuous(breaks=c(0,0.01,0.03,0.05,0.10)) +
  labs(x="MAF threshold", y="Mean FST") +
  theme_bw()

pdf("figures/Figure5.pdf", width=6, height=4)
print(p5)
dev.off()
cat("Figure 5 saved\n")

# ---- Figure 6: LD inflation ----
ld_infl <- data.frame(
  Platform=c("44K","HDRA"),
  Unpruned=c(0.609,0.512),
  Baseline=c(0.423,0.326)
)
ld_m <- melt(ld_infl, id.vars="Platform", variable.name="LD_pruning", value.name="FST")
ld_m$Platform <- factor(ld_m$Platform, levels=c("44K","HDRA"))

p6 <- ggplot(ld_m, aes(x=Platform, y=FST, fill=LD_pruning)) +
  geom_bar(stat="identity", position="dodge", alpha=0.85) +
  geom_text(aes(label=round(FST,3)), position=position_dodge(0.9), vjust=-0.3, size=3.5) +
  scale_fill_manual(values=c(Unpruned="#D73027", Baseline="#4575B4")) +
  labs(y="Mean FST") +
  annotate("text", x=0.7, y=0.55, label="1.44x", size=4, fontface="bold") +
  annotate("text", x=1.7, y=0.55, label="1.57x", size=4, fontface="bold") +
  theme_bw()

pdf("figures/Figure6.pdf", width=5, height=4)
print(p6)
dev.off()
cat("Figure 6 saved\n")

# ---- Figure 7: ADMIXTURE ----
q5 <- as.matrix(read.table("/home/yehia/sensitivity/rice_pruned3_K5.Q"))
ord5 <- order(apply(q5,1,which.max))
q2 <- as.matrix(read.table("/home/yehia/sensitivity/rice_pruned3_K2.Q"))
ord2 <- order(apply(q2,1,which.max))

pdf("figures/Figure7.pdf", width=10, height=4)
par(mfrow=c(1,2), mar=c(2,4,3,1))
barplot(t(q2[ord2,]), col=brewer.pal(3,"Set1")[1:2], border=NA, space=0,
        main="ADMIXTURE K = 2", ylab="Ancestry", xlab="Individuals")
barplot(t(q5[ord5,]), col=brewer.pal(6,"Set1"), border=NA, space=0,
        main="ADMIXTURE K = 5", ylab="Ancestry", xlab="Individuals")
dev.off()
cat("Figure 7 saved\n")

# ---- Figure 10: LD decay ----
tryCatch({
  ld_files_44 <- c("sens_D4_LDnone_ld.ld","sens_D4_LD08_ld.ld","sens_D4_LD05_ld.ld")
  ld_paths_44 <- file.path("/home/yehia/sensitivity", ld_files_44)
  ld_paths_44 <- ld_paths_44[file.exists(ld_paths_44)]
  
  cols_44 <- c("#2166AC","#67A9CF","#D1E5F0")
  
  pdf("figures/Figure10.pdf", width=10, height=5)
  par(mfrow=c(1,2), mar=c(4,4,3,1))
  
  # Left: 44K
  plot(NULL, NULL, xlim=c(1,5000), ylim=c(0,0.5), log="x",
       xlab="Distance (kb, log)", ylab=expression(Mean~r^2), main="44K LD decay (chr1)")
  abline(h=0.1, lty=2, col="gray")
  for(i in seq_along(ld_paths_44)) {
    d <- read.table(ld_paths_44[i], header=TRUE)
    d$dist <- abs(d$BP_B - d$BP_A)/1000
    d$bin <- cut(d$dist, c(0,10,50,100,200,500,1000,5000))
    agg <- aggregate(R2 ~ bin, d, mean)
    # Get bin midpoints dynamically
    bins <- strsplit(gsub("[^0-9,]","",as.character(agg$bin)), ",")
    agg$mid <- sapply(bins, function(x) mean(as.numeric(x)))
    lines(agg$mid, agg$R2, type="b", col=cols_44[i], pch=16, lwd=1.5)
  }
  legend("topright", legend=gsub("sens_D4_|_ld\\.ld","",basename(ld_paths_44)),
         fill=cols_44[seq_along(ld_paths_44)], cex=0.7, bty="n")
  
  # Right: HDRA
  ld_files_hd <- c("hdra_L_none_ld.ld","hdra_L_08_ld.ld","hdra_L_02_ld.ld")
  ld_paths_hd <- file.path("/home/yehia/hdra_rdp1", ld_files_hd)
  ld_paths_hd <- ld_paths_hd[file.exists(ld_paths_hd)]
  cols_hd <- c("#B2182B","#EF8A62","#FDDBC7")
  
  plot(NULL, NULL, xlim=c(1,5000), ylim=c(0,0.5), log="x",
       xlab="Distance (kb, log)", ylab=expression(Mean~r^2), main="HDRA LD decay (chr1)")
  abline(h=0.1, lty=2, col="gray")
  for(i in seq_along(ld_paths_hd)) {
    d <- read.table(ld_paths_hd[i], header=TRUE)
    d$dist <- abs(d$BP_B - d$BP_A)/1000
    d$bin <- cut(d$dist, c(0,10,50,100,200,500,1000,5000))
    agg <- aggregate(R2 ~ bin, d, mean)
    bins <- strsplit(gsub("[^0-9,]","",as.character(agg$bin)), ",")
    agg$mid <- sapply(bins, function(x) mean(as.numeric(x)))
    lines(agg$mid, agg$R2, type="b", col=cols_hd[i], pch=16, lwd=1.5)
  }
  legend("topright", legend=gsub("hdra_|_ld\\.ld","",basename(ld_paths_hd)),
         fill=cols_hd[seq_along(ld_paths_hd)], cex=0.7, bty="n")
  
  dev.off()
  cat("Figure 10 saved\n")
}, error=function(e) cat("Figure 10 not generated:", conditionMessage(e), "\n"))

cat("\nDone. Figures in figures/\n")
