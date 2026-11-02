# ======================================================================
# XENIUM TMA CASE STUDY · PART 2: downstream statistics (R)
# OmicLearningHub · Module 3 · https://githubrudramani.github.io/OmicLearningHub/xenium.html#downstream
#
# Reads the labelled object saved by xenium_tma_case_study.R and runs:
#   A. Cell composition statistics, ANIMAL as the replicate (never cells)
#      A1. treatment effect over time, adjusted for slide
#      A2. lesion progression within the same animal (paired)
#   B. Spatial analysis
#      B1. which cell types sit next to each other (neighbourhood enrichment)
#      B2. immune infiltration and the epithelium-stroma interface (distances)
#      B3. tissue layering: cell types by distance from the basal layer
#   C. Gene-level changes
#      C1. pseudobulk DE per cell type (with the same designs as A1 / A2)
#      C2. pathway tests (e.g. p53, DNA repair, inflammation) on panel genes
#
# Every test uses one value per animal (or per animal x condition). Column
# names and group names are settings, so the same code fits other designs.
# ======================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
  library(Matrix); library(dbscan); library(limma)
})
set.seed(1)

# ---------------------------- SETTINGS --------------------------------
OBJ_RDS <- "xenium_tma_results/runs/<run>/xenium_tma_seurat.rds"   # from the case-study script
OUT_DIR <- "xenium_tma_downstream"

# metadata columns in the object
COL_ANIMAL <- "mouse_id"; COL_GROUP <- "condition"; COL_TIME <- "stratum"
COL_SLIDE  <- "slide";    COL_CORE  <- "sample_id"; COL_CT   <- "celltype"

# A1 / C1: treatment effect over time (different animals per group)
TRT_REF <- "control"; TRT_TEST <- "treated"
TIME_LEVELS <- c("week4", "week8", "week16")     # in order; NULL = use the data as they are

# A2 / C1: progression within animals (paired). First level = reference.
PROG_LEVELS <- c("early_lesion", "advanced_lesion")
PROG_MERGE  <- NULL     # e.g. list(advanced_lesion = c("stage2", "stage3")) to pool stages first

# B: spatial (cell-type names from your labels)
EPI_TYPES    <- c("Basal_KC", "Differentiated_KC", "HairFollicle_KC")   # epithelium / tumour
STROMA_TYPES <- c("Fibroblast", "Endothelial", "Pericyte_SMC")
IMMUNE_TYPES <- c("T_cell", "Macrophage", "Macrophage_DC", "Dendritic_cell", "Mast_cell", "Neutrophil")
BASAL_TYPES  <- c("Basal_KC")                     # reference layer for B3
INTERFACE_DIST <- 30                               # um: "at the interface" = within this of both sides
LAYER_BINS <- c(0, 25, 50, 100, 200, Inf)          # um from the basal layer
NHOOD_K <- 10; NPERM <- 200                        # neighbours per cell; label shuffles per core

# C: gene level
PATHWAYS <- c("HALLMARK_P53_PATHWAY", "HALLMARK_DNA_REPAIR", "HALLMARK_UV_RESPONSE_UP",
              "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_G2M_CHECKPOINT")
SPECIES <- "Mus musculus"
MIN_ANIMALS <- 2; MIN_CELLS_PB <- 20
# ----------------------------------------------------------------------

# ========================== FUNCTIONS =================================
# One sample per animal x group: counts of each level of 'var' per sample,
# plus a sample table (time and slide taken from that animal's cores).
per_animal_counts <- function(md, var) {
  md$.s <- paste(md[[COL_ANIMAL]], md[[COL_GROUP]], sep = "|")
  cnt <- as.matrix(as.data.frame.matrix(table(md$.s, md[[var]])))
  mode1 <- function(x) names(sort(table(x), decreasing = TRUE))[1]
  samp <- md %>% group_by(.s) %>%
    summarise(animal = first(.data[[COL_ANIMAL]]), group = first(.data[[COL_GROUP]]),
              time = mode1(.data[[COL_TIME]]), slide = mode1(.data[[COL_SLIDE]]),
              n_slides = n_distinct(.data[[COL_SLIDE]]), n_cells = n(), .groups = "drop") %>%
    as.data.frame()
  if (any(samp$n_slides > 1)) message("Note: ", sum(samp$n_slides > 1), " sample(s) have cores on two slides; ",
                                      "their most common slide is used as the slide covariate.")
  rownames(samp) <- samp$.s
  list(cnt = cnt[samp$.s, , drop = FALSE], samp = samp)
}

# logit-transformed proportions + limma (speckle / propeller idea)
logit_limma <- function(cnt, design) {
  p <- (cnt + 0.5) / rowSums(cnt + 0.5)
  eBayes(lmFit(t(log(p / (1 - p))), design), robust = TRUE)
}
mean_props <- function(cnt, by) {
  p <- cnt / rowSums(cnt)
  t(sapply(split(seq_len(nrow(p)), by), function(i) round(colMeans(p[i, , drop = FALSE]), 4)))
}

# A1: treatment effect overall (adjusted for slide and time), per time point,
# and whether it changes over time (interaction).
composition_treatment <- function(cnt, samp) {
  samp$group <- factor(samp$group, levels = c(TRT_REF, TRT_TEST))
  samp$time  <- factor(samp$time, levels = if (is.null(TIME_LEVELS)) unique(samp$time) else TIME_LEVELS)
  samp$slide <- factor(samp$slide)
  print(table(time = samp$time, group = samp$group))          # animals per cell of the design
  adj <- if (nlevels(samp$slide) > 1) "slide + " else ""
  # overall effect
  d1 <- model.matrix(as.formula(paste("~", adj, "time + group")), samp)
  overall <- topTable(logit_limma(cnt, d1), coef = paste0("group", TRT_TEST), number = Inf)
  # per time point: one coefficient per time x group, contrasts within each time
  samp$tg <- factor(make.names(paste(samp$time, samp$group, sep = "_")))
  d2 <- model.matrix(as.formula(paste("~ 0 + tg", if (nzchar(adj)) "+ slide" else "")), samp)
  colnames(d2) <- sub("^tg", "", colnames(d2))
  cons <- c()
  for (t in levels(samp$time)) {
    a <- make.names(paste(t, TRT_TEST, sep = "_")); b <- make.names(paste(t, TRT_REF, sep = "_"))
    if (all(c(a, b) %in% colnames(d2))) cons[paste0(TRT_TEST, "_vs_", TRT_REF, "_at_", t)] <- paste(a, "-", b)
  }
  cmat <- makeContrasts(contrasts = cons, levels = d2); colnames(cmat) <- names(cons)
  fit2 <- eBayes(contrasts.fit(lmFit(t(qlogis((cnt + 0.5) / rowSums(cnt + 0.5))), d2), cmat), robust = TRUE)
  per_time <- lapply(names(cons), function(k) cbind(contrast = k, celltype = rownames(fit2),
                                                     topTable(fit2, coef = k, number = Inf, sort.by = "none")))
  # interaction: does the treatment effect differ between time points? (F-test)
  d3 <- model.matrix(as.formula(paste("~", adj, "time * group")), samp)
  int_cols <- grep(":", colnames(d3), value = TRUE)
  interaction <- if (length(int_cols)) topTable(logit_limma(cnt, d3), coef = int_cols, number = Inf) else NULL
  list(overall = cbind(celltype = rownames(overall), overall),
       per_time = do.call(rbind, per_time), interaction = interaction,
       mean_props = mean_props(cnt, paste(samp$time, samp$group)))
}

# A2: paired progression: same animal at several stages, animal as a block.
composition_paired <- function(cnt, samp, levels) {
  samp$group <- factor(samp$group, levels = levels)
  both <- names(which(tapply(samp$group, samp$animal, function(g) levels[1] %in% g && length(unique(g)) > 1)))
  keep <- samp$animal %in% both
  message("Paired animals: ", length(both))
  samp <- droplevels(samp[keep, ]); cnt <- cnt[keep, , drop = FALSE]
  d <- model.matrix(~ animal + group, samp)
  fit <- logit_limma(cnt, d)
  out <- lapply(levels[-1], function(l) if (paste0("group", l) %in% colnames(d))
    cbind(contrast = paste0(l, "_vs_", levels[1]), celltype = rownames(fit),
          topTable(fit, coef = paste0("group", l), number = Inf, sort.by = "none")))
  list(results = do.call(rbind, out), n_pairs = length(both), mean_props = mean_props(cnt, samp$group))
}

# B1: neighbourhood enrichment within one core. z > 0: type B is found next
# to type A more often than if labels were shuffled within the core.
nhood_enrichment_core <- function(xy, ct, k = 10, nperm = 200) {
  ct <- factor(ct); n <- length(ct); k <- min(k, n - 1)
  nn <- dbscan::kNN(xy, k = k)$id
  obs_mat <- function(lab) {
    m <- table(factor(rep(lab, k), levels = levels(ct)), factor(lab[as.vector(nn)], levels = levels(ct)))
    as.matrix(m)
  }
  obs  <- obs_mat(ct)
  perm <- replicate(nperm, obs_mat(sample(ct)))
  mu <- apply(perm, 1:2, mean); sdv <- apply(perm, 1:2, sd)
  z <- (obs - mu) / pmax(sdv, 1e-9)
  z[!is.finite(z)] <- NA
  z
}

# B2: distance from each 'from' cell to the nearest 'to' cell (within one core)
dist_to_nearest <- function(xy_from, xy_to) {
  if (!nrow(xy_from) || !nrow(xy_to)) return(rep(NA_real_, nrow(xy_from)))
  dbscan::kNN(xy_to, k = 1, query = xy_from)$dist[, 1]
}

# Per-core spatial summaries: infiltration, interface composition, layering
spatial_core_summary <- function(d) {
  xy <- as.matrix(d[, c("x_centroid", "y_centroid")]); ct <- d[[COL_CT]]
  epi <- ct %in% EPI_TYPES; str <- ct %in% STROMA_TYPES; imm <- ct %in% IMMUNE_TYPES
  d_epi <- dist_to_nearest(xy, xy[epi, , drop = FALSE])
  d_str <- dist_to_nearest(xy, xy[str, , drop = FALSE])
  interface <- d_epi <= INTERFACE_DIST & d_str <= INTERFACE_DIST
  bas <- ct %in% BASAL_TYPES
  d_bas <- dist_to_nearest(xy, xy[bas, , drop = FALSE])
  list(cells = data.frame(cell = rownames(d), celltype = ct, d_epi = d_epi, d_stroma = d_str,
                          d_basal = d_bas, interface = interface),
       summary = data.frame(
         n_immune = sum(imm),
         immune_within_epi = mean(d_epi[imm] <= INTERFACE_DIST, na.rm = TRUE),   # infiltration
         median_immune_to_epi = median(d_epi[imm], na.rm = TRUE),
         interface_cells = sum(interface, na.rm = TRUE),
         immune_frac_at_interface = mean(imm[which(interface)])))
}

# C1: pseudobulk counts per cell type x animal x group
pseudobulk <- function(counts, md) {
  id <- paste(md[[COL_CT]], md[[COL_ANIMAL]], md[[COL_GROUP]], sep = "||")
  f <- factor(id); mm <- sparse.model.matrix(~ 0 + f); colnames(mm) <- levels(f)
  pb <- as.matrix(counts[, rownames(md)] %*% mm)
  info <- data.frame(id = colnames(pb), do.call(rbind, strsplit(colnames(pb), "||", fixed = TRUE)),
                     n_cells = as.vector(table(f)[colnames(pb)]))
  names(info)[2:4] <- c("celltype", "animal", "group")
  s <- md %>% group_by(.a = .data[[COL_ANIMAL]], .g = .data[[COL_GROUP]]) %>%
    summarise(time = first(.data[[COL_TIME]]), slide = names(sort(table(.data[[COL_SLIDE]]), decreasing = TRUE))[1],
              .groups = "drop")
  k <- match(paste(info$animal, info$group), paste(s$.a, s$.g))
  info$time <- s$time[k]; info$slide <- s$slide[k]; rownames(info) <- info$id
  list(pb = pb, info = info)
}
# DESeq2 for one cell type with a given design; returns NULL + reason if under-powered
deseq_celltype <- function(pb, info, ct, design, contrast, levels) {
  ii <- info[info$celltype == ct & info$n_cells >= MIN_CELLS_PB & info$group %in% levels, ]
  n <- sapply(levels, function(l) length(unique(ii$animal[ii$group == l])))
  if (any(n < MIN_ANIMALS)) return(list(res = NULL, reason = paste(levels, n, sep = "=", collapse = ", ")))
  ii$group <- factor(ii$group, levels = levels); ii$animal <- factor(ii$animal)
  ii$time <- factor(ii$time); ii$slide <- factor(ii$slide)
  for (v in c("slide", "time")) if (nlevels(ii[[v]]) < 2) design <- update(design, as.formula(paste("~ . -", v)))
  dds <- DESeq2::DESeqDataSetFromMatrix(round(pb[, ii$id]), ii, design)
  dds <- dds[rowSums(DESeq2::counts(dds) >= 5) >= 2, ]
  dds <- tryCatch(DESeq2::DESeq(dds, quiet = TRUE), error = function(e) e)
  if (inherits(dds, "error")) return(list(res = NULL, reason = conditionMessage(dds)))
  r <- as.data.frame(DESeq2::results(dds, contrast = contrast)); r$gene <- rownames(r)
  list(res = r[order(r$pvalue), ], reason = NA, design = deparse(design), n = n)
}
# C2: pathway test on one DE result: fgsea on the Wald statistic, gene sets
# cut to the panel genes that were tested (never the whole genome)
pathway_test <- function(res, sets, min_size = 5) {
  st <- res$stat; names(st) <- res$gene; st <- sort(st[!is.na(st)], decreasing = TRUE)
  sets <- lapply(sets, intersect, names(st))
  sets <- sets[lengths(sets) >= min_size]
  if (!length(sets)) return(NULL)
  r <- fgsea::fgsea(sets, st, minSize = min_size, maxSize = 500)
  r$size_on_panel <- lengths(sets)[r$pathway]
  r[order(r$pval), ]
}
# Per-cell pathway score: mean of the z-scored genes of the set (panel genes only)
pathway_score <- function(lognorm, genes) {
  g <- intersect(genes, rownames(lognorm)); if (length(g) < 3) return(rep(NA_real_, ncol(lognorm)))
  x <- as.matrix(lognorm[g, , drop = FALSE]); x <- (x - rowMeans(x)) / pmax(apply(x, 1, sd), 1e-9)
  colMeans(x)
}
# =======================================================================

# ============================== RUN ===================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
obj <- readRDS(OBJ_RDS)
md  <- obj@meta.data
md  <- md[!is.na(md[[COL_ANIMAL]]) & !md[[COL_CT]] %in% c("Unassigned", "Ambiguous", "Low_quality_undefined"), ]
if (!is.null(PROG_MERGE)) for (nm in names(PROG_MERGE)) md[[COL_GROUP]][md[[COL_GROUP]] %in% PROG_MERGE[[nm]]] <- nm

# ---- A1. treatment effect over time, adjusted for slide ----------------
# Proportions add up to 1: if one cell type expands, every other type's share
# falls even when its cell numbers do not change. Read a significant drop in
# a large population next to the rise that caused it, and confirm key results
# with ratios (e.g. T cells per epithelial cell) or cells per mm2 of tissue.
a <- md[md[[COL_GROUP]] %in% c(TRT_REF, TRT_TEST), ]
if (nrow(a)) {
  pa <- per_animal_counts(a, COL_CT)
  A1 <- composition_treatment(pa$cnt, pa$samp)
  write.csv(A1$overall,  file.path(OUT_DIR, "A1_composition_treatment_overall.csv"), row.names = FALSE)
  write.csv(A1$per_time, file.path(OUT_DIR, "A1_composition_treatment_per_time.csv"), row.names = FALSE)
  if (!is.null(A1$interaction)) write.csv(cbind(celltype = rownames(A1$interaction), A1$interaction),
                                          file.path(OUT_DIR, "A1_composition_treatment_x_time.csv"), row.names = FALSE)
  pp <- as.data.frame(pa$cnt / rowSums(pa$cnt)); pp$time <- pa$samp$time; pp$group <- pa$samp$group
  ggsave(file.path(OUT_DIR, "A1_composition_over_time.pdf"), width = 14, height = 9,
         ggplot(pivot_longer(pp, -c(time, group), names_to = "celltype", values_to = "prop"),
                aes(factor(time, levels = if (is.null(TIME_LEVELS)) unique(time) else TIME_LEVELS), prop, colour = group)) +
           geom_point(position = position_dodge(0.5)) +
           stat_summary(aes(group = group), fun = mean, geom = "line", position = position_dodge(0.5)) +
           facet_wrap(~ celltype, scales = "free_y") + labs(x = NULL, y = "Proportion per animal") + theme_bw())
}

# ---- A2. progression within the same animal (paired) -------------------
p <- md[md[[COL_GROUP]] %in% PROG_LEVELS, ]
if (nrow(p)) {
  pp <- per_animal_counts(p, COL_CT)
  A2 <- composition_paired(pp$cnt, pp$samp, PROG_LEVELS)
  if (!is.null(A2$results)) write.csv(A2$results, file.path(OUT_DIR, "A2_composition_progression_paired.csv"), row.names = FALSE)
}

# ---- B1. neighbourhood enrichment, per core -> per animal ----------------
z_all <- list()
for (core in unique(md[[COL_CORE]])) {
  d <- md[md[[COL_CORE]] == core, ]; if (nrow(d) < 50) next
  z <- nhood_enrichment_core(as.matrix(d[, c("x_centroid", "y_centroid")]), d[[COL_CT]], NHOOD_K, NPERM)
  zz <- as.data.frame(as.table(z)); names(zz) <- c("from", "to", "z")
  zz$core <- core; zz$animal <- d[[COL_ANIMAL]][1]; zz$group <- d[[COL_GROUP]][1]; zz$time <- d[[COL_TIME]][1]
  z_all[[core]] <- zz
}
z_all <- do.call(rbind, z_all)
z_animal <- z_all %>% group_by(animal, group, time, from, to) %>% summarise(z = mean(z, na.rm = TRUE), .groups = "drop")
write.csv(z_animal, file.path(OUT_DIR, "B1_neighbourhood_enrichment_per_animal.csv"), row.names = FALSE)
ggsave(file.path(OUT_DIR, "B1_neighbourhood_enrichment_by_group.pdf"), width = 14, height = 6,
       ggplot(z_animal %>% group_by(group, from, to) %>% summarise(z = mean(z, na.rm = TRUE), .groups = "drop"),
              aes(to, from, fill = pmax(pmin(z, 10), -10))) + geom_tile() + facet_wrap(~ group) +
         scale_fill_gradient2(low = "#2563EB", high = "#C2410C", name = "z") + RotatedAxis() + theme_bw() +
         labs(x = "neighbour", y = "cell", title = "Neighbourhood enrichment (mean z per animal)"))

# ---- B2 / B3. infiltration, interface, layering ---------------------------
sp <- lapply(split(md, md[[COL_CORE]]), spatial_core_summary)
cells_sp <- do.call(rbind, lapply(sp, `[[`, "cells"))
core_sum <- do.call(rbind, lapply(names(sp), function(k) cbind(core = k, sp[[k]]$summary)))
cm <- md %>% distinct(core = .data[[COL_CORE]], animal = .data[[COL_ANIMAL]], group = .data[[COL_GROUP]],
                      time = .data[[COL_TIME]], slide = .data[[COL_SLIDE]])
core_sum <- merge(core_sum, cm, by = "core")
animal_sum <- core_sum %>% group_by(animal, group, time, slide) %>%
  summarise(across(c(immune_within_epi, median_immune_to_epi, immune_frac_at_interface),
                   ~ weighted.mean(.x, n_immune, na.rm = TRUE)), n_immune = sum(n_immune), .groups = "drop")
write.csv(animal_sum, file.path(OUT_DIR, "B2_infiltration_interface_per_animal.csv"), row.names = FALSE)
# simple test per metric: animal-level values, same adjustments as A1
tr <- animal_sum[animal_sum$group %in% c(TRT_REF, TRT_TEST), ]
if (length(unique(tr$group)) == 2) {
  tr$group <- factor(tr$group, levels = c(TRT_REF, TRT_TEST))
  f <- paste("~", if (length(unique(tr$slide)) > 1) "slide +" else "", if (length(unique(tr$time)) > 1) "time +" else "", "group")
  b2 <- do.call(rbind, lapply(c("immune_within_epi", "median_immune_to_epi", "immune_frac_at_interface"), function(v) {
    s <- summary(lm(as.formula(paste(v, f)), data = tr))$coefficients
    data.frame(metric = v, effect = s[paste0("group", TRT_TEST), 1], p = s[paste0("group", TRT_TEST), 4])
  }))
  write.csv(b2, file.path(OUT_DIR, "B2_infiltration_tests.csv"), row.names = FALSE)
}
cells_sp$animal <- md[cells_sp$cell, COL_ANIMAL]; cells_sp$group <- md[cells_sp$cell, COL_GROUP]
layers <- cells_sp %>% filter(!celltype %in% BASAL_TYPES, !is.na(d_basal)) %>%
  mutate(bin = cut(d_basal, LAYER_BINS, include.lowest = TRUE)) %>%
  count(group, animal, bin, celltype) %>% group_by(group, animal, bin) %>% mutate(frac = n / sum(n)) %>% ungroup()
write.csv(layers, file.path(OUT_DIR, "B3_layers_from_basal.csv"), row.names = FALSE)
ggsave(file.path(OUT_DIR, "B3_layers_from_basal.pdf"), width = 12, height = 5,
       ggplot(layers %>% group_by(group, bin, celltype) %>% summarise(frac = mean(frac), .groups = "drop"),
              aes(bin, frac, fill = celltype)) + geom_col() + facet_wrap(~ group) +
         labs(x = "distance from the basal layer (um)", y = "mean fraction per animal") + theme_bw() + RotatedAxis())

# ---- C1. pseudobulk DE per cell type ---------------------------------------
pbl <- pseudobulk(LayerData(obj, "counts", assay = "Xenium"), md)
sets <- tryCatch({
  m <- if ("collection" %in% names(formals(msigdbr::msigdbr))) msigdbr::msigdbr(species = SPECIES, collection = "H")
       else msigdbr::msigdbr(species = SPECIES, category = "H")
  split(m$gene_symbol, m$gs_name)[PATHWAYS]
}, error = function(e) { message("msigdbr failed (", conditionMessage(e), ") - pathway steps skipped"); NULL })
de_runs <- list(
  treatment = list(levels = c(TRT_REF, TRT_TEST), design = ~ slide + time + group),
  progression = list(levels = PROG_LEVELS[1:2], design = ~ animal + group))   # add more pairs if needed
summ <- list()
for (nm in names(de_runs)) {
  run <- de_runs[[nm]]; dir.create(file.path(OUT_DIR, paste0("C_", nm)), showWarnings = FALSE)
  for (ct in sort(unique(pbl$info$celltype))) {
    r <- deseq_celltype(pbl$pb, pbl$info, ct, run$design, c("group", run$levels[2], run$levels[1]), run$levels)
    if (is.null(r$res)) { summ[[paste(nm, ct)]] <- data.frame(analysis = nm, celltype = ct, tested = FALSE,
                                                              note = paste("animals:", r$reason)); next }
    write.csv(r$res, file.path(OUT_DIR, paste0("C_", nm), paste0("DE_", make.names(ct), ".csv")), row.names = FALSE)
    pw <- if (!is.null(sets)) pathway_test(r$res, sets) else NULL
    if (!is.null(pw)) data.table::fwrite(pw, file.path(OUT_DIR, paste0("C_", nm), paste0("pathways_", make.names(ct), ".tsv")),
                                         sep = "\t", sep2 = c("", ",", ""))
    summ[[paste(nm, ct)]] <- data.frame(analysis = nm, celltype = ct, tested = TRUE, note = r$design,
                                        de_genes_padj05 = sum(r$res$padj < 0.05, na.rm = TRUE))
  }
}
write.csv(bind_rows(summ), file.path(OUT_DIR, "C_DE_summary.csv"), row.names = FALSE)

# ---- C2. pathway scores per cell, averaged per animal and cell type -------
if (!is.null(sets)) {
  ln <- LayerData(obj, "data", assay = "Xenium")[, rownames(md)]
  sc <- sapply(sets, function(g) pathway_score(ln, g))
  scd <- cbind(md[, c(COL_ANIMAL, COL_GROUP, COL_TIME, COL_CT)], sc)
  names(scd)[1:4] <- c("animal", "group", "time", "celltype")
  sca <- scd %>% group_by(animal, group, time, celltype) %>%
    summarise(across(all_of(names(sets)), ~ mean(.x, na.rm = TRUE)), n = n(), .groups = "drop")
  write.csv(sca, file.path(OUT_DIR, "C2_pathway_scores_per_animal.csv"), row.names = FALSE)
  ggsave(file.path(OUT_DIR, "C2_pathway_scores.pdf"), width = 16, height = 10,
         ggplot(pivot_longer(sca, all_of(names(sets)), names_to = "pathway", values_to = "score"),
                aes(group, score, colour = group)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15) +
           facet_grid(pathway ~ celltype, scales = "free_y") + theme_bw() + RotatedAxis() +
           theme(strip.text.y = element_text(angle = 0, size = 6)))
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Done. Results in ", normalizePath(OUT_DIR))
