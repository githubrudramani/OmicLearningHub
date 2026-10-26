# ======================================================================
# OmicLearningHub · Module 1 · complete, reproducible run
# Bulk RNA-seq with limma-voom + TMM on GSE63310 (Law et al. 2016),
# followed by MSigDB GSEA / ORA (scripts/rnaseq_gsea_msigdb.R).
#
#   Rscript run_rnaseq_limma.R [output_dir]      (default: output)
#
# Needs internet on the first run (GEO data, MSigDB gene sets).
# Every result the module shows is written to output_dir, and
# checks_vs_paper.csv compares key numbers with the published workflow.
# ======================================================================
t0 <- Sys.time()
suppressPackageStartupMessages({
  library(limma); library(edgeR); library(org.Mm.eg.db); library(ggplot2)
  library(pheatmap); library(data.table)
})
args    <- commandArgs(trailingOnly = TRUE)
OUT     <- normalizePath(if (length(args)) args[1] else "output", mustWork = FALSE)
here    <- tryCatch(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))),
                    error = function(e) getwd())
GSEA_R  <- normalizePath(file.path(here, "..", "scripts", "rnaseq_gsea_msigdb.R"), mustWork = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "plots"), showWarnings = FALSE)
options(timeout = 1800)
set.seed(1)
checks <- list()
check <- function(name, observed, expected, tol = 0) {
  ok <- if (is.numeric(expected)) all(abs(observed - expected) <= tol) else identical(observed, expected)
  checks[[length(checks) + 1]] <<- data.frame(check = name,
    expected = paste(expected, collapse = " "), observed = paste(observed, collapse = " "),
    status = if (ok) "match" else "DIFFERS")
  message(sprintf("[%s] %s", if (ok) "match" else "DIFFERS", name))
}
png_plot <- function(file, expr, w = 1600, h = 1200) {
  png(file.path(OUT, "plots", file), width = w, height = h, res = 200); on.exit(dev.off()); expr
}

# 03 · download GSE63310 and build a DGEList ----------------------------
dl <- file.path(OUT, "data"); dir.create(dl, showWarnings = FALSE)
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt",    "GSM1545540_JMS8-3.txt",   "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt",    "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
if (!all(file.exists(file.path(dl, files)))) {
  tarf <- file.path(dl, "GSE63310_RAW.tar")
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file",
                destfile = tarf, mode = "wb")
  untar(tarf, exdir = dl)
  for (f in file.path(dl, paste0(files, ".gz"))) R.utils::gunzip(f, overwrite = TRUE)
}
x <- readDGE(file.path(dl, files), columns = c(1, 3))
colnames(x) <- substring(colnames(x), 12)
check("genes in raw count table", nrow(x), 27179)

# 04 · samples, batch, annotation ---------------------------------------
group <- factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))
lane  <- factor(rep(c("L004", "L006", "L008"), c(3, 4, 2)))
x$samples$group <- group; x$samples$lane <- lane
genes <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = rownames(x),
                          columns = c("SYMBOL", "GENENAME"), keytype = "ENTREZID"))
genes <- genes[!duplicated(genes$ENTREZID), ]
x$genes <- genes

# 05 · filter ------------------------------------------------------------
lcpm_raw <- cpm(x, log = TRUE)
keep <- filterByExpr(x, group = group)
x <- x[keep, , keep.lib.sizes = FALSE]
check("genes kept by filterByExpr", nrow(x), 16624, tol = 166)      # within 1%
png_plot("05_density_before_after_filtering.png", {
  par(mfrow = c(1, 2)); plotDensities(lcpm_raw, legend = FALSE, main = "A. Raw")
  plotDensities(cpm(x, log = TRUE), legend = FALSE, main = "B. Filtered") }, w = 2000, h = 900)

# 06 · TMM ---------------------------------------------------------------
x <- calcNormFactors(x, method = "TMM")
check("TMM normalisation factors", round(x$samples$norm.factors, 3),
      c(0.894, 1.025, 1.046, 1.046, 1.016, 0.922, 0.996, 1.086, 0.984), tol = 0.002)
lcpm <- cpm(x, log = TRUE)
write.csv(x$samples, file.path(OUT, "06_samples_libsize_normfactors.csv"))

# 07 · MDS ----------------------------------------------------------------
col.group <- c("#C2410C", "#2563EB", "#16865A")[as.integer(group)]
png_plot("07_MDS.png", { par(mfrow = c(1, 2))
  plotMDS(lcpm, labels = group, col = col.group, main = "A. Cell type")
  plotMDS(lcpm, labels = lane, dim.plot = c(3, 4), main = "B. Lane (dims 3 & 4)") }, w = 2000, h = 1000)

# 08 · design & contrasts -----------------------------------------------
design <- model.matrix(~ 0 + group + lane)
colnames(design) <- gsub("group", "", colnames(design))
contr.matrix <- makeContrasts(BasalvsLP = Basal - LP, BasalvsML = Basal - ML, LPvsML = LP - ML,
                              levels = colnames(design))

# 09 · voom → lmFit → eBayes ----------------------------------------------
png_plot("09_voom_mean_variance.png", v <- voom(x, design, plot = TRUE))
vfit <- contrasts.fit(lmFit(v, design), contrasts = contr.matrix)
efit <- eBayes(vfit)
png_plot("09_SA_plot.png", plotSA(efit, main = "Final model: mean-variance trend"))
de <- summary(decideTests(efit))
write.csv(unclass(de), file.path(OUT, "09_decideTests_eBayes.csv"))
check("eBayes DE counts (Down/NotSig/Up x 3 contrasts)", as.vector(unclass(de)),
      c(4648, 7113, 4863, 4927, 7026, 4671, 3135, 10972, 2517), tol = 100)

# 10 · treat, tables, plots ----------------------------------------------
tfit <- treat(vfit, lfc = 1)
dt <- decideTests(tfit)
write.csv(unclass(summary(dt)), file.path(OUT, "10_decideTests_treat_lfc1.csv"))
check("treat |log2FC|>1 DE counts", as.vector(unclass(summary(dt))),
      c(1632, 12976, 2016, 1777, 12790, 2057, 224, 16210, 190), tol = 60)
for (k in colnames(contr.matrix))
  write.csv(topTreat(tfit, coef = k, n = Inf), file.path(OUT, paste0("10_DE_", k, "_treat.csv")), row.names = FALSE)
write.fit(tfit, dt, file = file.path(OUT, "10_all_contrasts_fit.txt"))
png_plot("10_MD_BasalvsLP.png",
         plotMD(tfit, column = 1, status = dt[, 1], main = colnames(tfit)[1], xlim = c(-8, 13)))
tt <- topTable(efit, coef = "BasalvsLP", n = Inf)
tt$status <- ifelse(tt$adj.P.Val < 0.05 & tt$logFC > 1, "Up",
             ifelse(tt$adj.P.Val < 0.05 & tt$logFC < -1, "Down", "NS"))
ggsave(file.path(OUT, "plots", "10_volcano_BasalvsLP.png"),
       ggplot(tt, aes(logFC, -log10(P.Value), colour = status)) + geom_point(size = 0.5) +
         scale_colour_manual(values = c(Up = "#C2410C", Down = "#2563EB", NS = "grey70")) +
         theme_minimal(), width = 7, height = 5, dpi = 200)
basal.vs.lp <- topTreat(tfit, coef = 1, n = Inf)
top <- rownames(basal.vs.lp)[1:100]
png(file.path(OUT, "plots", "10_heatmap_top100_BasalvsLP.png"), width = 1400, height = 1800, res = 200)
pheatmap(lcpm[top, ], scale = "row", show_rownames = FALSE,
         annotation_col = data.frame(group = group, row.names = colnames(lcpm)))
dev.off()

# 11 · camera with MSigDB Hallmark ---------------------------------------
suppressPackageStartupMessages(library(msigdbr))
h <- if ("collection" %in% names(formals(msigdbr))) msigdbr(species = "Mus musculus", collection = "H") else
       msigdbr(species = "Mus musculus", category = "H")
idcol <- intersect(c("ncbi_gene", "entrez_gene"), names(h))[1]
idx <- ids2indices(split(as.character(h[[idcol]]), h$gs_name), id = rownames(v))
cam <- camera(v, idx, design, contrast = contr.matrix[, "BasalvsLP"])
write.csv(cam, file.path(OUT, "11_camera_Hallmark_BasalvsLP.csv"))
png_plot("11_barcode_EMT_BasalvsLP.png",
         barcodeplot(efit$t[, "BasalvsLP"], index = idx[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
                     main = "Basal vs LP: EMT hallmark"))

# 12–13 · GSEA (all genes ranked by t) + ORA / restricted GSEA -----------
if (file.exists(GSEA_R)) {
  owd <- setwd(OUT); source(GSEA_R, local = new.env()); setwd(owd)
} else message("GSEA script not found at ", GSEA_R, " – skipping steps 12–13")

# record --------------------------------------------------------------------
chk <- do.call(rbind, checks)
write.csv(chk, file.path(OUT, "checks_vs_paper.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
writeLines(c(sprintf("run_date: %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
             sprintf("runtime_minutes: %s", mins),
             sprintf("R: %s", R.version.string),
             sprintf("platform: %s", R.version$platform),
             sprintf("checks_matching: %d of %d", sum(chk$status == "match"), nrow(chk))),
           file.path(OUT, "run_summary.txt"))
print(chk, row.names = FALSE)
message("Done in ", mins, " min. Results in ", OUT)
