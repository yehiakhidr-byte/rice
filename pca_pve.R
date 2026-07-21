library(adegenet)
obj <- read.PLINK("rice_pruned3.raw", parallel=FALSE)
gt <- tab(obj)
for(i in 1:ncol(gt)){gt[is.na(gt[,i]), i] <- mean(gt[,i], na.rm=TRUE)}
pca <- prcomp(gt, scale=FALSE, center=TRUE)
pve <- summary(pca)$importance[2,1:4]*100
cat(sprintf("PC1: %.1f%%\nPC2: %.1f%%\nPC3: %.1f%%\nPC4: %.1f%%\n", pve[1], pve[2], pve[3], pve[4]))
