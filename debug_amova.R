SENS_DIR <- "/home/yehia/sensitivity"
OUT_DIR <- "/mnt/d/project genomic analysis"

library(vegan)

# Use rice_pruned3 as test case
sc <- "rice_pruned3"
mdist_file <- file.path(SENS_DIR, paste0(sc, "_dist.mdist"))
mdist_id_file <- file.path(SENS_DIR, paste0(sc, "_dist.mdist.id"))

# Read ID file
ids <- read.table(mdist_id_file, header=FALSE, stringsAsFactors=FALSE)[,2]
cat("IDs:", length(ids), "\n")

# Read distance matrix - try as-is
m <- as.matrix(read.table(mdist_file, header=FALSE, sep="\t"))
cat("Matrix dim:", dim(m), "\n")
cat("Match IDs:", nrow(m) == length(ids), "\n")

# Give row/col names
rownames(m) <- colnames(m) <- ids

# Create population labels (first 5 as group1, rest as group2)
pop_f <- factor(rep(c("A","B"), c(5, nrow(m)-5)))

cat("Pop levels:", levels(pop_f), "\n")
cat("Running adonis2...\n")

set.seed(999)
ad <- adonis2(as.dist(m) ~ pop_f, permutations=99)
cat("Result:\n")
print(ad)
cat("R2:", ad["pop_f", "R2"], "\n")
cat("F:", ad["pop_f", "F"], "\n")
cat("Pr:", ad["pop_f", "Pr(>F)"], "\n")
