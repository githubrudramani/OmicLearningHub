# Install every package used by Module 1 into an existing R (>= 4.5) without conda.
# Run once:  Rscript install.R
cran <- c("BiocManager", "R.utils", "ggplot2", "pheatmap", "msigdbr", "data.table")
bioc <- c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "fgsea", "clusterProfiler")
new <- setdiff(cran, rownames(installed.packages()))
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE, ask = FALSE)
for (p in c(cran, bioc)) cat(sprintf("%-16s %s\n", p, as.character(packageVersion(p))))
