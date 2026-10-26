# ======================================================================
# GSEA WITH MSigDB AFTER limma-voom
# OmicLearningHub · Module 1 · https://githubrudramani.github.io/OmicLearningHub/rnaseq.html
#
# Run after the Module 1 pipeline (needs `efit` from eBayes, whose genes
# table has a SYMBOL column). Two complementary analyses per collection:
#
#   A) Ranked GSEA, ALL genes ranked by limma's moderated t-statistic
#      (fgsea). Uses every gene, no cut-off; finds coordinated shifts.
#   B) SIGNIFICANT genes only (adj.P.Val < FDR_CUTOFF):
#      B1) over-representation analysis (ORA) of up and down genes
#          separately, against all tested genes as the universe
#          (clusterProfiler::enricher) – the standard test for a gene list;
#      B2) GSEA restricted to the significant genes – tests only whether a
#          set's DE genes lean up or down. Exploratory; read the note below.
# ======================================================================

suppressPackageStartupMessages({
  library(limma); library(fgsea); library(clusterProfiler); library(msigdbr)
  library(data.table); library(ggplot2)
})
set.seed(42)

# ---------------- SETTINGS --------------------------------------------
CONTRAST   <- "BasalvsLP"          # a column of contr.matrix
SPECIES    <- "Mus musculus"       # "Homo sapiens" for human data
FDR_CUTOFF <- 0.05                 # significant genes for part B
MIN_SIZE   <- 15; MAX_SIZE <- 500  # gene-set size limits (after matching)
OUT_DIR    <- "gsea_results"
COLLECTIONS <- list(               # name = c(collection, subcollection)
  Hallmark = c("H",  NA),
  Reactome = c("C2", "CP:REACTOME"),
  GO_BP    = c("C5", "GO:BP"))
# ----------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE)

# 1. GENE SETS FROM MSigDB ----------------------------------------------
# msigdbr >= 10 uses collection/subcollection; older versions use
# category/subcategory. Recent versions download the gene sets on first use,
# so the first run needs internet access.
get_msigdb <- function(species, coll, subcoll) {
  sub <- if (is.na(subcoll)) NULL else subcoll
  m <- if ("collection" %in% names(formals(msigdbr))) {
    msigdbr(species = species, collection = coll, subcollection = sub)
  } else {
    msigdbr(species = species, category = coll, subcategory = sub)
  }
  unique(as.data.frame(m)[, c("gs_name", "gene_symbol")])
}
sets <- lapply(COLLECTIONS, function(cl) get_msigdb(SPECIES, cl[1], cl[2]))
print(sapply(sets, function(s) length(unique(s$gs_name))))   # sets per collection

# 2. RANK ALL GENES BY THE MODERATED t ----------------------------------
tt <- topTable(efit, coef = CONTRAST, number = Inf, sort.by = "none")
tt$SYMBOL <- if ("SYMBOL" %in% colnames(tt)) as.character(tt$SYMBOL) else rownames(tt)
tt <- tt[!is.na(tt$SYMBOL) & tt$SYMBOL != "", ]
tt <- tt[order(-abs(tt$t)), ]
tt <- tt[!duplicated(tt$SYMBOL), ]                    # one row per symbol (strongest t)
ranks    <- sort(setNames(tt$t, tt$SYMBOL), decreasing = TRUE)
universe <- names(ranks)                              # every tested gene
sig_up   <- tt$SYMBOL[tt$adj.P.Val < FDR_CUTOFF & tt$logFC > 0]
sig_down <- tt$SYMBOL[tt$adj.P.Val < FDR_CUTOFF & tt$logFC < 0]
message(length(ranks), " ranked genes; significant (adj.P < ", FDR_CUTOFF, "): ",
        length(sig_up), " up, ", length(sig_down), " down")

save_fgsea <- function(res, file) {                   # leadingEdge is a list column
  res <- res[order(res$padj, -abs(res$NES)), ]
  fwrite(res, file, sep = "\t", sep2 = c("", ",", ""))
  res
}
nes_plot <- function(res, title, n = 10) {
  d <- res[!is.na(res$padj) & res$padj < 0.05, ]
  d <- rbind(head(d[order(-d$NES), ], n), head(d[order(d$NES), ], n))
  d <- d[!duplicated(d$pathway), ]                  # unique() fails on the leadingEdge list column
  if (nrow(d) == 0) return(NULL)
  d$pathway <- factor(d$pathway, levels = d$pathway[order(d$NES)])
  ggplot(d, aes(NES, pathway, fill = NES > 0)) + geom_col() +
    scale_fill_manual(values = c(`TRUE` = "#C2410C", `FALSE` = "#2563EB"),
                      labels = c(`TRUE` = "Up", `FALSE` = "Down"), name = NULL) +
    labs(title = title, x = "Normalised enrichment score (NES)", y = NULL) +
    theme_minimal(base_size = 10)
}

results <- list()
for (nm in names(sets)) {
  pw  <- split(sets[[nm]]$gene_symbol, sets[[nm]]$gs_name)
  t2g <- sets[[nm]]

  # A) RANKED GSEA, ALL GENES ---------------------------------------------
  gsea_all <- fgsea(pw, ranks, minSize = MIN_SIZE, maxSize = MAX_SIZE)
  gsea_all <- save_fgsea(gsea_all, file.path(OUT_DIR, paste0(CONTRAST, "_", nm, "_GSEA_allgenes.tsv")))
  p <- nes_plot(gsea_all, paste(nm, "· GSEA, all genes ranked by t ·", CONTRAST))
  if (!is.null(p)) ggsave(file.path(OUT_DIR, paste0(CONTRAST, "_", nm, "_GSEA_allgenes.pdf")),
                          p, width = 9, height = 6)
  top <- gsea_all$pathway[which.min(gsea_all$padj)]
  if (length(top) == 1)
    ggsave(file.path(OUT_DIR, paste0(CONTRAST, "_", nm, "_top_enrichment_plot.pdf")),
           plotEnrichment(pw[[top]], ranks) + labs(title = top), width = 7, height = 4)

  # B1) ORA OF SIGNIFICANT GENES (up and down separately) -----------------
  ora <- function(genes, dir) {
    if (length(genes) < 10) return(NULL)
    e <- enricher(genes, TERM2GENE = t2g, universe = universe,
                  minGSSize = MIN_SIZE, maxGSSize = MAX_SIZE,
                  pvalueCutoff = 1, qvalueCutoff = 1)       # keep all; filter on p.adjust later
    if (is.null(e)) return(NULL)
    d <- as.data.frame(e); if (nrow(d) == 0) return(NULL)
    d$direction <- dir; d
  }
  ora_res <- rbind(ora(sig_up, "up"), ora(sig_down, "down"))
  if (!is.null(ora_res)) {
    ora_res <- ora_res[order(ora_res$p.adjust), ]
    fwrite(ora_res, file.path(OUT_DIR, paste0(CONTRAST, "_", nm, "_ORA_sig_genes.tsv")), sep = "\t")
  }

  # B2) GSEA RESTRICTED TO SIGNIFICANT GENES (exploratory) ----------------
  # Every gene here is DE, so this does not ask "is the set changed?" but
  # "do the set's DE genes lean up or down?". Sets must still contain
  # >= 10 significant genes, so many sets drop out.
  ranks_sig <- ranks[names(ranks) %in% c(sig_up, sig_down)]
  gsea_sig <- fgsea(pw, ranks_sig, minSize = 10, maxSize = MAX_SIZE)
  gsea_sig <- save_fgsea(gsea_sig, file.path(OUT_DIR, paste0(CONTRAST, "_", nm, "_GSEA_siggenes.tsv")))

  results[[nm]] <- list(gsea_all = gsea_all, ora = ora_res, gsea_sig = gsea_sig)
  message(sprintf("%-9s GSEA all genes: %3d sets padj<0.05 | ORA up/down: %3d/%3d | GSEA sig-only: %3d",
                  nm, sum(gsea_all$padj < 0.05, na.rm = TRUE),
                  if (is.null(ora_res)) 0L else sum(ora_res$p.adjust < 0.05 & ora_res$direction == "up"),
                  if (is.null(ora_res)) 0L else sum(ora_res$p.adjust < 0.05 & ora_res$direction == "down"),
                  sum(gsea_sig$padj < 0.05, na.rm = TRUE)))
}

# 3. COMPARE THE APPROACHES (Hallmark) ------------------------------------
h <- results$Hallmark
cmp <- merge(h$gsea_all[, c("pathway", "NES", "padj")],
             h$gsea_sig[, c("pathway", "NES", "padj")],
             by = "pathway", all = TRUE, suffixes = c(".all", ".sig"))
if (!is.null(h$ora)) {
  o <- h$ora[, c("ID", "direction", "p.adjust")]
  o <- reshape(o, idvar = "ID", timevar = "direction", direction = "wide")
  cmp <- merge(cmp, o, by.x = "pathway", by.y = "ID", all = TRUE)
}
cmp <- cmp[order(cmp$padj.all), ]
fwrite(cmp, file.path(OUT_DIR, paste0(CONTRAST, "_Hallmark_method_comparison.tsv")), sep = "\t")
print(as.data.frame(head(cmp, 15)), row.names = FALSE)

saveRDS(results, file.path(OUT_DIR, paste0(CONTRAST, "_gsea_results.rds")))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Done. Results in ", normalizePath(OUT_DIR))
