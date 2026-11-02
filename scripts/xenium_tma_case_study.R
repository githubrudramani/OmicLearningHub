# ======================================================================
# XENIUM TMA CASE STUDY · multi-slide tissue microarray pipeline (R)
# OmicLearningHub · Module 3 · https://githubrudramani.github.io/OmicLearningHub/xenium.html#case-study
#
# Example design: two Xenium slides, each a tissue microarray (TMA) of mouse
# tissue cores. Every animal contributes one or more cores; every group being
# compared is spread over both slides. Marker lists use mouse skin on the
# Xenium Mouse 5K panel as the example: replace them for your tissue.
#
# Backbone: 10x Xenium 5K R tutorial (median-scaled normalisation, sketch
# 25k cells per slide -> cluster -> project to all cells) and the Seurat v5
# spatial vignette (image plots, RCTD, niches).
#
# Safety rules built in
#   * Replicate = animal. Cores from one animal are never separate replicates;
#     paired analyses (same animal in both conditions) use complete pairs only.
#   * Statistics (composition test, pseudobulk DE) only run when every kept
#     core has a verified mouse_id + condition in curated/core_map.csv AND the
#     cluster labels are marked reviewed in curated/cluster_labels.csv.
#     Until then the run produces exploratory outputs only.
#   * core_map.csv stores each core's centre; cluster_labels.csv stores a
#     fingerprint of the clustering. If cores move or the clustering changes,
#     the old file is set aside (_OLD_<date>.csv), never silently reused.
#   * If a compared group sits on only one slide, group and slide effects are
#     inseparable: the script says so and the result is exploratory.
#   * Every run writes to its own folder runs/<RUN_NAME>/; nothing is deleted.
#     Hand-edited files live in curated/ and are read by every run.
#
# Run top to bottom, or section by section in RStudio / Jupyter.
# ======================================================================

library(Seurat); library(arrow); library(dbscan)
library(ggplot2); library(dplyr); library(tidyr); library(Matrix)
options(future.globals.maxSize = 50 * 1024^3)  # limit for parallel workers, not a RAM allocation
set.seed(1)

# ---------------------------- SETTINGS --------------------------------
BASE_DIR     <- path.expand("path/to/project")          # folder holding the Xenium output bundles
RESULTS_ROOT <- path.expand("xenium_tma_results")
#   curated/       your hand-edited inputs: tma_layout_*.csv, tma_cases.csv,
#                  core_map.csv, cluster_labels.csv
#   runs/<name>/   one folder per run, never overwritten by later runs
RUN_NAME   <- format(Sys.time(), "run_%Y%m%d_%H%M%S")    # or a fixed name, e.g. "run_qc20"
SLIDE_DIRS <- c(TMA1 = "xenium_output_slide1",          # Xenium output folders inside BASE_DIR
                TMA2 = "xenium_output_slide2")

# Study design (columns of curated/core_map.csv)
#   mouse_id  = animal; every core of that animal shares it, also across
#               tissues and slides
#   condition = what you compare (e.g. control / treated, or tissue type)
#   stratum   = optional second factor (e.g. time point); "none" if unused
COMPARE_REF  <- "control"     # names must match core_map.csv exactly
COMPARE_TEST <- "treated"
DE_SUBSET_STRATUM <- NA       # compare within one stratum, e.g. "week8" (NA = all)
PAIRED   <- FALSE             # TRUE if the same animal gives both conditions
                              # (e.g. lesion vs. adjacent normal tissue)
MIN_MICE <- 2                 # per group per cell type (>= 3 recommended)
MIN_CELLS_PB <- 20            # cells per pseudobulk sample

# Core detection (microns)
CORE_EPS <- 75; CORE_MINPTS <- 10
CORE_MIN_CELLS  <- 300        # smaller groups are treated as debris (check the core map!)
CORE_MERGE_DIST <- 800        # merge fragments whose centres are closer than this;
                              # keep it below the centre spacing of neighbouring cores
ROW_GAP <- 800                # y jump between TMA rows, for provisional R01C01 names

# TMA map (optional, recommended). Put these in curated/:
#   tma_layout_<slide>.csv  one row per map position:
#                           core_number, block_id ("Empty" for empty spots), map_row, map_col
#   tma_cases.csv           one row per block: slide, block_id, mouse_id, condition, stratum,
#                           plus any per-animal columns you like (sex, age, ...)
# With both, cores are matched to the map automatically and core_map.csv is pre-filled.
TMA_LAYOUT_FILES <- c(TMA1 = "tma_layout_TMA1.csv", TMA2 = "tma_layout_TMA2.csv")
CASES_FILE       <- "tma_cases.csv"
TMA_ORIENTATION  <- c(TMA1 = NA, TMA2 = NA)  # NA = automatic; 1-8 forces one (see cores/*_orientations.csv)

# Start over: TRUE deletes curated/core_map.csv and cluster_labels.csv at the
# start of the run (a copy is kept in the run folder). Use it once, then FALSE.
RESET_CURATED <- FALSE

# QC (starting values from the 10x tutorial; check qc/ outputs before trusting)
QC_MIN_COUNTS <- 20; QC_MAX_COUNT_Q <- 0.98   # upper cut-off computed PER SLIDE (NA = off)
QC_MIN_GENES  <- 15
QC_MAX_CTRL_RATIO <- 0.05    # (control probe + codeword counts) / gene counts
QC_MIN_AREA <- 10; QC_MAX_AREA_Q <- 0.995     # per slide

# Annotation
ANNOT_MIN_SCORE  <- 0        # best module score must exceed this, else "Unassigned"
ANNOT_MIN_MARGIN <- 0.05     # best - second best below this -> "Ambiguous"

# Niches
NICHE_K <- 30; N_NICHES <- 8
NICHE_MAX_DIST <- 100        # ignore neighbours further than this (um); neighbours stay within a core

GENES_SHOW <- c("Trp63", "Dsg1a", "Krt79", "Col1a1", "Ptprc", "Pecam1", "Mki67")  # example genes to map

# Optional steps
RUN_ENRICHR <- FALSE         # 8a: Enrichr cell-type lookup per cluster (internet)
RUN_RCTD    <- FALSE         # 8b: reference labels (spacexr); build the reference with
REF_RDS     <- "path/to/reference.rds"   #     scripts/build_skin_rctd_reference.R
REF_CT_COL  <- "celltype"
RCTD_ALL_CELLS <- FALSE      # FALSE = sketch cells only
RUN_ZOOM    <- TRUE; ZOOM_N <- 4
RUN_BANKSY  <- FALSE
# ----------------------------------------------------------------------

# ======================= HELPER FUNCTIONS ==============================
# (kept together so they can be read, reused and tested on their own)

# DBSCAN finds dense groups of cells; fragments whose centres are closer than
# merge_dist are merged; groups < min_cells are dropped. Returns per-cell group
# ids (0 = not in a core) and provisional R..C.. names by TMA row / column.
detect_cores <- function(xy, eps, minpts, merge_dist, min_cells, row_gap) {
  cl <- dbscan::dbscan(xy, eps = eps, minPts = minpts)$cluster
  if (!any(cl > 0)) return(NULL)
  cen <- aggregate(xy[cl > 0, , drop = FALSE], list(cl = cl[cl > 0]), median)
  if (nrow(cen) > 1) {
    grp <- cutree(hclust(dist(cen[, 2:3]), method = "single"), h = merge_dist)
    cl[cl > 0] <- grp[match(cl[cl > 0], cen$cl)]
  }
  sz <- table(cl[cl > 0])
  cl[!cl %in% as.integer(names(sz)[sz >= min_cells])] <- 0
  if (!any(cl > 0)) return(NULL)
  cc <- aggregate(xy[cl > 0, , drop = FALSE], list(cl = cl[cl > 0]), median)
  names(cc)[2:3] <- c("x", "y")
  cc <- cc[order(cc$y), ]
  cc$row <- cumsum(c(TRUE, diff(cc$y) > row_gap))
  cc <- cc[order(cc$row, cc$x), ]
  cc$col <- ave(cc$x, cc$row, FUN = seq_along)
  cc$core_id <- sprintf("R%02dC%02d", cc$row, cc$col)
  list(cell_core = cc$core_id[match(cl, cc$cl)], cores = cc)
}

# Split sorted coordinates at the k-1 largest gaps -> group index (1..k).
# 'sep' = smallest chosen gap / largest remaining gap (higher = clearer grid).
.gap_groups <- function(v, k) {
  if (k == 1) return(list(grp = rep(1L, length(v)), sep = Inf))
  if (k > length(v)) return(NULL)
  o <- order(v); g <- diff(v[o])
  cut_at <- sort(order(g, decreasing = TRUE)[seq_len(k - 1)])
  grp <- integer(length(v)); grp[o] <- findInterval(seq_along(v), cut_at + 1) + 1L
  rest <- g[-cut_at]
  list(grp = grp, sep = min(g[cut_at]) / max(c(rest, 1e-9)))
}

# Match detected core centres (x, y) to the TMA map. The grid is inferred from
# the centres; the 8 rotations / mirror images are tried and the one whose
# "Empty" map positions fall on missing cores wins. Slides are often scanned
# rotated or flipped relative to the paper map, so never assume orientation.
match_tma_grid <- function(x, y, lay) {
  is_empty <- grepl("^empty$", trimws(lay$block_id), ignore.case = TRUE)
  Rm <- max(lay$map_row); Cm <- max(lay$map_col)
  L <- matrix(NA_integer_, Rm, Cm); L[cbind(lay$map_row, lay$map_col)] <- as.integer(lay$core_number)
  best <- NULL
  for (dm in list(c(Rm, Cm), c(Cm, Rm))) {
    gy <- .gap_groups(y, dm[1]); gx <- .gap_groups(x, dm[2])
    if (is.null(gy) || is.null(gx) || anyDuplicated(paste(gy$grp, gx$grp))) next
    sc <- min(gy$sep, gx$sep)
    if (is.null(best) || sc > best$sep) best <- list(R = dm[1], C = dm[2], row = gy$grp, col = gx$grp, sep = sc)
  }
  if (is.null(best)) return(NULL)
  empty_nums <- as.integer(lay$core_number[is_empty]); full_nums <- as.integer(lay$core_number[!is_empty])
  res <- list(); nums <- list()
  for (tr in c(FALSE, TRUE)) for (fr in c(FALSE, TRUE)) for (fc in c(FALSE, TRUE)) {
    a <- if (tr) best$col else best$row; b <- if (tr) best$row else best$col
    Rl <- if (tr) best$C else best$R;    Cl <- if (tr) best$R else best$C
    if (fr) a <- Rl + 1 - a
    if (fc) b <- Cl + 1 - b
    ok_shape <- nrow(L) == Rl && ncol(L) == Cl
    num <- if (ok_shape) L[cbind(a, b)] else NA
    mism <- if (ok_shape) sum(num %in% empty_nums) + sum(!(full_nums %in% num)) else Inf
    res[[length(res) + 1]] <- data.frame(transpose = tr, flip_rows = fr, flip_cols = fc, mismatches = mism)
    nums[[length(nums) + 1]] <- num
  }
  res <- do.call(rbind, res)
  list(grid = best, table = res, nums = nums, best = which(res$mismatches == min(res$mismatches)))
}

# Validate core_map.csv against the current cores. Returns TRUE only if every
# row is usable; prints what is wrong otherwise.
check_core_map <- function(core_map, cur_cores, move_tol, paired) {
  ok <- TRUE
  if (anyDuplicated(core_map$sample_id)) stop("Duplicated sample_id in core_map.csv")
  ex <- toupper(trimws(as.character(core_map$exclude)))
  if (!all(ex %in% c("TRUE", "FALSE"))) {
    message("WARNING: 'exclude' must be TRUE or FALSE for every core; bad rows: ",
            paste(core_map$sample_id[!ex %in% c("TRUE", "FALSE")], collapse = ", "))
    ok <- FALSE
  }
  chk <- merge(cur_cores, core_map[, c("sample_id", "x_center", "y_center")], by = "sample_id",
               all.x = TRUE, suffixes = c("", "_map"))
  moved <- chk$sample_id[!is.na(chk$x_center_map) &
             sqrt((chk$x_center - chk$x_center_map)^2 + (chk$y_center - chk$y_center_map)^2) > move_tol]
  if (length(moved)) {
    message("WARNING: these core IDs now sit somewhere else than in core_map.csv (core detection changed?): ",
            paste(moved, collapse = ", "), ". Rebuild core_map.csv before using its metadata.")
    ok <- FALSE
  }
  keep <- ex == "FALSE"
  mc <- unique(core_map[keep & !is.na(core_map$mouse_id), c("mouse_id", "condition", "stratum")])
  if (!paired && anyDuplicated(mc$mouse_id)) {
    message("WARNING: a mouse_id has more than one condition/stratum but PAIRED = FALSE: ",
            paste(unique(mc$mouse_id[duplicated(mc$mouse_id)]), collapse = ", "))
    ok <- FALSE
  }
  ok
}

# Fingerprint of a clustering: which sketch cell is in which cluster.
clustering_fingerprint <- function(cells, clusters) {
  f <- tempfile(); on.exit(unlink(f))
  writeLines(paste(sort(paste(cells, clusters)), collapse = ";"), f)
  unname(tools::md5sum(f))
}

# Neighbourhood composition per cell, computed WITHIN each core and only up to
# max_dist, so no cell borrows neighbours from another core.
niche_profiles <- function(md, ct_levels, k, max_dist) {
  out <- list(); n_nb <- setNames(integer(nrow(md)), rownames(md))
  for (sid in unique(md$sample_id)) {
    d <- md[md$sample_id == sid, ]; n <- nrow(d); if (n < 2) next
    kk <- min(k, n - 1)
    nn <- dbscan::kNN(as.matrix(d[, c("x_centroid", "y_centroid")]), k = kk)
    id <- nn$id; id[nn$dist > max_dist] <- NA
    ii <- rep(seq_len(n), kk); jj <- as.vector(id); ok <- !is.na(jj)
    ct <- factor(d$celltype, levels = ct_levels)
    m  <- as.matrix(sparseMatrix(i = ii[ok], j = as.integer(ct[jj[ok]]), x = 1,
                                 dims = c(n, length(ct_levels)), dimnames = list(rownames(d), ct_levels)))
    nv <- rowSums(m); n_nb[rownames(d)] <- nv
    out[[sid]] <- (m / pmax(nv, 1))[nv > 0, , drop = FALSE]
  }
  list(mat = if (length(out)) do.call(rbind, out) else NULL, n_neighbours = n_nb)
}

# Animals present in every group (for paired designs).
complete_pairs <- function(df, lv) {
  tab <- table(df$mouse_id, factor(df$grp, levels = lv))
  df[df$mouse_id %in% rownames(tab)[rowSums(tab > 0) == length(lv)], , drop = FALSE]
}

# Composition test: logit-transformed proportions + limma (the speckle /
# propeller approach); one sample = animal x condition; animal blocked if paired.
composition_test <- function(md2, lv, paired, min_mice) {
  samp <- md2 %>% distinct(rep_id, mouse_id, grp) %>% as.data.frame()
  if (paired) samp <- complete_pairs(samp, lv)
  n_m <- sapply(lv, function(g) length(unique(samp$mouse_id[samp$grp == g])))
  if (any(n_m < min_mice)) { message("Composition test skipped: animals per group ",
                                     paste(lv, n_m, sep = "=", collapse = ", ")); return(NULL) }
  cnt <- as.matrix(as.data.frame.matrix(table(md2$rep_id, md2$celltype)))[samp$rep_id, , drop = FALSE]
  props <- (cnt + 0.5) / rowSums(cnt + 0.5)
  Y <- t(log(props / (1 - props)))
  samp$mouse_id <- factor(samp$mouse_id)
  des <- if (paired) model.matrix(~ mouse_id + grp, samp) else model.matrix(~ grp, samp)
  fit <- limma::eBayes(limma::lmFit(Y, des), robust = TRUE)
  tt  <- limma::topTable(fit, coef = paste0("grp", lv[2]), number = Inf, sort.by = "P")
  tt$celltype <- rownames(tt)
  for (g in lv) tt[[paste0("mean_prop_", g)]] <- round(colMeans(props[samp$grp == g, tt$celltype, drop = FALSE]), 4)
  tt
}

# Pseudobulk DE per cell type: counts summed per cell type x animal x condition,
# DESeq2 with ~ grp (or ~ mouse_id + grp when paired).
pseudobulk_de <- function(counts, md2, lv, paired, min_mice, min_cells) {
  md2$pb_id <- paste(md2$celltype, md2$mouse_id, md2$condition, sep = "||")
  fac <- factor(md2$pb_id)
  mm  <- sparse.model.matrix(~ 0 + fac); colnames(mm) <- levels(fac)
  pb  <- as.matrix(counts[, rownames(md2)] %*% mm)
  info <- md2 %>% group_by(pb_id, celltype, mouse_id, condition, grp) %>%
    summarise(n_cells = n(), .groups = "drop") %>% as.data.frame()
  if (anyDuplicated(info$pb_id)) stop("An animal has several groups within one condition - check core_map.csv")
  res <- list(); skipped <- list(); summ <- list()
  for (ct in sort(unique(info$celltype))) {
    if (ct %in% c("Unassigned", "Ambiguous")) {
      skipped[[ct]] <- data.frame(celltype = ct, reason = "uncertain label (mixture) - not tested"); next }
    ii <- info[info$celltype == ct & info$n_cells >= min_cells, ]
    if (paired) ii <- complete_pairs(ii, lv)
    n_m <- sapply(lv, function(g) length(unique(ii$mouse_id[ii$grp == g])))
    if (any(n_m < min_mice)) {
      skipped[[ct]] <- data.frame(celltype = ct, reason = paste0(if (paired) "complete pairs" else "animals",
        " with >= ", min_cells, " cells: ", paste(lv, n_m, sep = "=", collapse = ", "))); next }
    ii$mouse_id <- droplevels(factor(ii$mouse_id)); rownames(ii) <- ii$pb_id
    design <- if (paired) ~ mouse_id + grp else ~ grp
    dds <- DESeq2::DESeqDataSetFromMatrix(round(pb[, ii$pb_id]), colData = ii, design = design)
    dds <- dds[rowSums(DESeq2::counts(dds) >= 5) >= 2, ]
    dds <- tryCatch(DESeq2::DESeq(dds, quiet = TRUE), error = function(e) e)
    if (inherits(dds, "error")) {
      skipped[[ct]] <- data.frame(celltype = ct, reason = paste("DESeq2 error:", conditionMessage(dds))); next }
    r <- as.data.frame(DESeq2::results(dds, contrast = c("grp", lv[2], lv[1])))
    r$gene <- rownames(r)
    res[[ct]] <- r[order(r$padj), c("gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
    summ[[ct]] <- data.frame(celltype = ct, design = deparse(design), t(setNames(n_m, paste0("n_mice_", lv))),
                             genes_tested = nrow(r), padj_lt_0.05 = sum(r$padj < 0.05, na.rm = TRUE))
  }
  list(results = res, info = info, summary = if (length(summ)) do.call(rbind, summ),
       skipped = if (length(skipped)) do.call(rbind, skipped))
}
# =======================================================================

dirs <- setNames(file.path(BASE_DIR, SLIDE_DIRS), names(SLIDE_DIRS))
stopifnot(all(dir.exists(dirs)))
CURATED_DIR <- file.path(RESULTS_ROOT, "curated")
OUT_DIR     <- file.path(RESULTS_ROOT, "runs", RUN_NAME)
if (dir.exists(OUT_DIR) && length(list.files(OUT_DIR)))
  stop("Run folder already exists and is not empty: ", OUT_DIR, " - choose another RUN_NAME.")
for (d in c("qc", "cores", "clusters", "spatial", "niches", "de", "explorer"))
  dir.create(file.path(OUT_DIR, d), showWarnings = FALSE, recursive = TRUE)
dir.create(CURATED_DIR, showWarnings = FALSE, recursive = TRUE)
writeLines(OUT_DIR, file.path(RESULTS_ROOT, "latest_run.txt"))
settings <- mget(ls(pattern = "^[A-Z][A-Z0-9_]+$"))
writeLines(c(paste("Run started:", format(Sys.time())),
             paste(names(settings), sapply(settings, function(x) paste(deparse(x), collapse = "")), sep = " = ")),
           file.path(OUT_DIR, "run_settings.txt"))

# 0. OPTIONAL RESET OF CURATED FILES -------------------------------------
if (RESET_CURATED) for (f in file.path(CURATED_DIR, c("core_map.csv", "cluster_labels.csv"))) if (file.exists(f)) {
  file.copy(f, file.path(OUT_DIR, paste0("deleted_", basename(f)))); file.remove(f)
  message("Deleted ", f, " (copy kept in the run folder). Set RESET_CURATED back to FALSE.")
}

# 1. LOAD SLIDES ---------------------------------------------------------
# No cell outlines (they can overflow R's protection stack on a full slide)
# and no transcript positions; step 9b reads outlines per core instead.
cell_cols <- c("x_centroid", "y_centroid", "cell_area", "nucleus_area", "nucleus_count",
               "segmentation_method", "control_probe_counts", "control_codeword_counts",
               "genomic_control_counts")
objs <- list()
for (s in names(dirs)) {
  o  <- LoadXenium(data.dir = dirs[[s]], fov = s, assay = "Xenium", molecule.coordinates = FALSE)
  cm <- as.data.frame(arrow::read_parquet(file.path(dirs[[s]], "cells.parquet")))
  rownames(cm) <- as.character(cm$cell_id)
  o  <- AddMetaData(o, cm[colnames(o), intersect(cell_cols, colnames(cm)), drop = FALSE])
  o$slide <- s
  objs[[s]] <- o
}
rm(o, cm); gc()

# 2. MERGE ---------------------------------------------------------------
obj <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))   # "Key 'Xenium_' taken" warning is harmless
rm(objs); gc()
DefaultAssay(obj) <- "Xenium"
obj[["Xenium"]] <- JoinLayers(obj[["Xenium"]])
obj$ctrl_ratio <- rowSums(cbind(obj$control_probe_counts, obj$control_codeword_counts),
                          na.rm = TRUE) / pmax(obj$nCount_Xenium, 1)

# 3. DETECT TMA CORES AND MATCH THEM TO THE MAP --------------------------
obj$core_id <- NA_character_
for (s in names(dirs)) {
  cells <- colnames(obj)[obj$slide == s]
  xy <- as.matrix(obj@meta.data[cells, c("x_centroid", "y_centroid")])
  dc <- detect_cores(xy, CORE_EPS, CORE_MINPTS, CORE_MERGE_DIST, CORE_MIN_CELLS, ROW_GAP)
  if (is.null(dc)) { message("WARNING: no cores found on ", s, " - check CORE_EPS / CORE_MIN_CELLS"); next }
  obj$core_id[match(cells, colnames(obj))] <- dc$cell_core
  d <- obj@meta.data[cells, ]
  lab <- d %>% filter(!is.na(core_id)) %>% group_by(core_id) %>%
    summarise(x = median(x_centroid), y = median(y_centroid), n = n(), .groups = "drop")
  ggsave(file.path(OUT_DIR, "cores", paste0(s, "_core_map.pdf")), width = 12, height = 12,
         ggplot(d[sample(nrow(d), min(nrow(d), 2e5)), ], aes(x_centroid, y_centroid, colour = is.na(core_id))) +
           geom_point(size = 0.05) + scale_y_reverse() + coord_fixed() +
           scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "red"), guide = "none") +
           geom_label(data = lab, aes(x, y, label = paste0(core_id, "\n", n)), inherit.aes = FALSE, size = 2) +
           theme_void() + ggtitle(paste(s, "- detected cores (red = not in a core)")))
}

obj$tma_core <- NA_integer_
for (s in names(dirs)) {
  lay_file <- file.path(CURATED_DIR, TMA_LAYOUT_FILES[[s]])
  if (is.na(TMA_LAYOUT_FILES[s]) || !file.exists(lay_file)) {
    message("No TMA layout for ", s, " - cores keep R..C.. names; fill core_map.csv by hand."); next }
  lay <- read.csv(lay_file, stringsAsFactors = FALSE)
  cc <- obj@meta.data[obj$slide == s & !is.na(obj$core_id), ] %>% group_by(core_id) %>%
    summarise(x = median(x_centroid), y = median(y_centroid), .groups = "drop")
  m <- match_tma_grid(cc$x, cc$y, lay)
  if (is.null(m)) { message("WARNING: ", s, ": could not fit the map grid to ", nrow(cc), " detected cores."); next }
  write.csv(cbind(orientation = seq_len(nrow(m$table)), m$table),
            file.path(OUT_DIR, "cores", paste0(s, "_orientations.csv")), row.names = FALSE)
  pick <- if (!is.na(TMA_ORIENTATION[s])) TMA_ORIENTATION[[s]] else m$best
  message(s, ": grid ", m$grid$R, " x ", m$grid$C, ", ", nrow(cc), " cores; orientation ",
          paste(pick, collapse = "/"), " with ", m$table$mismatches[pick[1]], " mismatches (0 = perfect)")
  if (length(pick) > 1) { message("WARNING: several orientations fit equally - set TMA_ORIENTATION[\"", s,
                                  "\"] after checking cores/", s, "_orientations.csv"); next }
  if (m$table$mismatches[pick] > 0) message("WARNING: ", s, ": ", m$table$mismatches[pick],
    " position(s) do not fit the map (lost core, or a fragment detected as a core).")
  num <- m$nums[[pick]]
  if (anyDuplicated(num)) { message("WARNING: ", s, ": two cores on one map position - not matched."); next }
  idx <- obj$slide == s & !is.na(obj$core_id)
  obj$tma_core[idx] <- num[match(obj$core_id[idx], cc$core_id)]
  obj$core_id[idx]  <- sprintf("P%02d", obj$tma_core[idx])            # stable map numbers
  cc$num <- num; cc$block <- lay$block_id[match(num, lay$core_number)]
  d <- obj@meta.data[obj$slide == s, ]
  ggsave(file.path(OUT_DIR, "cores", paste0(s, "_core_map_matched.pdf")), width = 14, height = 12,
         ggplot(d[sample(nrow(d), min(nrow(d), 2e5)), ], aes(x_centroid, y_centroid)) +
           geom_point(size = 0.05, colour = "grey50") + scale_y_reverse() + coord_fixed() +
           geom_label(data = cc, aes(x, y, label = paste0("#", num, "\n", block)), size = 2) +
           theme_void() + ggtitle(paste(s, "- cores matched to the TMA map (check against the H&E!)")))
}
obj$sample_id <- ifelse(is.na(obj$core_id), NA, paste(obj$slide, obj$core_id, sep = "_"))

# Core -> animal / condition. The template is pre-filled from the TMA map +
# cases table where cores were matched; everything else stays EMPTY on
# purpose. Each core's centre is stored so a later run can tell if an ID now
# points to a different physical core.
core_map_file <- file.path(CURATED_DIR, "core_map.csv")
map_cols <- c("slide", "core_id", "sample_id", "n_cells", "x_center", "y_center",
              "mouse_id", "condition", "stratum", "exclude")
cur_cores <- obj@meta.data %>% filter(!is.na(sample_id)) %>% group_by(slide, core_id, sample_id) %>%
  summarise(n_cells = n(), x_center = round(median(x_centroid)), y_center = round(median(y_centroid)),
            .groups = "drop")
if (file.exists(core_map_file) &&
    !all(map_cols %in% names(read.csv(core_map_file, nrows = 1, check.names = FALSE)))) {
  bak <- sub("\\.csv$", paste0("_OLD_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"), core_map_file)
  file.rename(core_map_file, bak); message("WARNING: core_map.csv was in an old format; saved as ", basename(bak))
}
if (!file.exists(core_map_file)) {
  tmpl <- cur_cores %>% mutate(mouse_id = NA_character_, condition = NA_character_, stratum = NA_character_,
                               exclude = FALSE, notes = "", block_id = NA_character_)
  tmpl$tma_core <- as.integer(ifelse(grepl("^P[0-9]+$", tmpl$core_id), sub("^P", "", tmpl$core_id), NA))
  cases <- if (file.exists(file.path(CURATED_DIR, CASES_FILE)))
    read.csv(file.path(CURATED_DIR, CASES_FILE), stringsAsFactors = FALSE) else NULL
  for (s in unique(tmpl$slide)) {
    lay_file <- file.path(CURATED_DIR, TMA_LAYOUT_FILES[[s]])
    r <- tmpl$slide == s & !is.na(tmpl$tma_core)
    if (!file.exists(lay_file) || !any(r)) next
    lay <- read.csv(lay_file, stringsAsFactors = FALSE)
    tmpl$block_id[r] <- lay$block_id[match(tmpl$tma_core[r], lay$core_number)]
    if (is.null(cases)) next
    ci <- match(paste(s, tmpl$block_id[r]), paste(cases$slide, cases$block_id))
    for (cl in setdiff(names(cases), c("slide", "block_id"))) {
      if (!cl %in% names(tmpl)) tmpl[[cl]] <- NA
      tmpl[[cl]][r] <- cases[[cl]][ci]
    }
    if (any(is.na(ci))) message("WARNING: blocks not in ", CASES_FILE, ": ",
                                paste(unique(tmpl$block_id[r][is.na(ci)]), collapse = ", "))
  }
  write.csv(tmpl, core_map_file, row.names = FALSE, na = "")
  message(">>> Wrote ", core_map_file, ". Check it against the core map PDFs, fill mouse_id / condition / ",
          "stratum for every core and set exclude = TRUE for damaged cores. Then re-run.")
}
core_map <- read.csv(core_map_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))
META_OK <- check_core_map(core_map, cur_cores, move_tol = CORE_EPS * 2, paired = PAIRED)
core_map$exclude <- toupper(trimws(as.character(core_map$exclude))) == "TRUE"
k <- match(obj$sample_id, core_map$sample_id)
obj$mouse_id  <- core_map$mouse_id[k]
obj$condition <- core_map$condition[k]
obj$stratum   <- core_map$stratum[k]
obj$exclude   <- core_map$exclude[k] %in% TRUE
reserved <- c(map_cols, "notes")
for (cl in setdiff(colnames(core_map), c(reserved, colnames(obj@meta.data))))   # per-animal extras
  obj[[cl]] <- core_map[[cl]][k]
write.csv(as.data.frame.matrix(table(condition = obj$condition, slide = obj$slide, useNA = "ifany")),
          file.path(OUT_DIR, "cores/condition_by_slide.csv"))      # are groups spread over slides?

# 4. QC -------------------------------------------------------------------
md0 <- obj@meta.data
max_count <- if (is.na(QC_MAX_COUNT_Q)) rep(Inf, nrow(md0)) else
  ave(md0$nCount_Xenium, md0$slide, FUN = function(x) quantile(x, QC_MAX_COUNT_Q))
max_area  <- ave(md0$cell_area, md0$slide, FUN = function(x) quantile(x, QC_MAX_AREA_Q, na.rm = TRUE))
fail <- data.frame(
  low_counts  = md0$nCount_Xenium < QC_MIN_COUNTS,
  high_counts = md0$nCount_Xenium > max_count,        # likely two cells segmented as one
  low_genes   = md0$nFeature_Xenium < QC_MIN_GENES,
  high_ctrl   = md0$ctrl_ratio > QC_MAX_CTRL_RATIO,
  small_cell  = md0$cell_area %in% NA | md0$cell_area < QC_MIN_AREA,
  large_cell  = !(md0$cell_area %in% NA) & md0$cell_area > max_area,
  not_in_core = is.na(md0$core_id),
  excluded    = md0$exclude, row.names = rownames(md0))
fail[is.na(fail)] <- FALSE
md0$qc_pass <- rowSums(fail) == 0
extra_qc <- intersect(c("genomic_control_counts", "nucleus_area", "nucleus_count"), colnames(md0))
qc_rep <- cbind(md0[, c("slide", "core_id", "nCount_Xenium", extra_qc)], fail, pass = md0$qc_pass) %>%
  filter(!is.na(core_id)) %>% group_by(slide, core_id) %>%
  summarise(n_before = n(), n_kept = sum(pass), pct_kept = round(100 * mean(pass), 1),
            across(low_counts:large_cell, sum), median_counts = median(nCount_Xenium),
            across(all_of(extra_qc), ~ round(median(.x, na.rm = TRUE), 2), .names = "median_{.col}"),
            .groups = "drop")
write.csv(qc_rep, file.path(OUT_DIR, "qc/qc_per_core.csv"), row.names = FALSE)   # why each core lost cells
pdf(file.path(OUT_DIR, "qc/qc_plots.pdf"), width = 10, height = 6)
sub <- md0[sample(nrow(md0), min(nrow(md0), 1e5)), ]
print(ggplot(sub, aes(nCount_Xenium, nFeature_Xenium, colour = qc_pass)) + geom_point(size = 0.1) +
        scale_x_log10() + scale_y_log10() + facet_wrap(~slide) + theme_bw())
print(ggplot(filter(md0, !is.na(core_id)), aes(core_id, nCount_Xenium)) + geom_boxplot(outlier.size = 0.1) +
        scale_y_log10() + facet_wrap(~slide, scales = "free_x", ncol = 1) + theme_bw() +
        theme(axis.text.x = element_text(angle = 90, size = 5)))
if ("nucleus_count" %in% extra_qc)               # 0 or 2+ nuclei = possible segmentation problem
  print(ggplot(sub, aes(factor(pmin(nucleus_count, 3)), nCount_Xenium)) + geom_violin() +
          scale_y_log10() + facet_wrap(~slide) + theme_bw() + xlab("nuclei per cell (3 = 3+)"))
for (s in names(dirs))                            # are failed cells localised (folds, edges)?
  print(ggplot(filter(sub, slide == s), aes(x_centroid, y_centroid, colour = qc_pass)) +
          geom_point(size = 0.05) + scale_y_reverse() + coord_fixed() + theme_void() + ggtitle(s))
dev.off()
obj <- subset(obj, cells = rownames(md0)[md0$qc_pass])
rm(md0, sub, fail); gc()

has_value <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))
kept <- obj@meta.data
incomplete <- unique(kept$sample_id[!(has_value(kept$mouse_id) & has_value(kept$condition) & has_value(kept$stratum))])
if (length(incomplete)) message("Kept cores with unfinished metadata: ", length(incomplete))
META_OK <- META_OK && nrow(kept) > 0 && !length(incomplete)
message("Metadata complete: ", META_OK); rm(kept)
obj[["Xenium"]] <- split(obj[["Xenium"]], f = obj$slide)

# 5. NORMALISE & SKETCH (10x 5K tutorial) ----------------------------------
obj <- NormalizeData(obj, scale.factor = median(obj$nCount_Xenium))   # median counts, as 10x
obj <- FindVariableFeatures(obj, nfeatures = 2000)
hv <- VariableFeatures(obj)          # keep genes seen in >= 10 cells of EVERY slide:
for (ly in Layers(obj[["Xenium"]], search = "counts"))   # an all-zero gene breaks LeverageScore
  hv <- hv[Matrix::rowSums(LayerData(obj, assay = "Xenium", layer = ly)[hv, , drop = FALSE] > 0) >= 10]
obj <- SketchData(obj, ncells = 25000, method = "LeverageScore", sketched.assay = "sketch", features = hv)

# 6. CLUSTER THE SKETCH -----------------------------------------------------
# No Harmony by default: check clusters/umap.pdf by slide first. If slides
# already overlap, Harmony only reshuffles clusters.
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, nfeatures = 2000) |> ScaleData() |> RunPCA(npcs = 30)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30) |> FindClusters(resolution = 0.5)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, return.model = TRUE, reduction.name = "umap")

# 7. PROJECT TO ALL CELLS ---------------------------------------------------
obj <- ProjectData(obj, assay = "Xenium", sketched.assay = "sketch", sketched.reduction = "pca",
                   full.reduction = "pca.full", dims = 1:30, umap.model = "umap",
                   refdata = list(cluster_full = "seurat_clusters"))
DefaultAssay(obj) <- "Xenium"
obj[["Xenium"]] <- JoinLayers(obj[["Xenium"]])
obj$cluster_full <- factor(obj$cluster_full, levels = sort(as.numeric(unique(obj$cluster_full))))
ggsave(file.path(OUT_DIR, "clusters/umap.pdf"), width = 14, height = 6,
       DimPlot(obj, reduction = "full.umap", group.by = c("cluster_full", "slide"), label = TRUE, raster = TRUE))
## seurat_clusters = sketch cells only (NA elsewhere); cluster_full = every cell

# 8. MARKERS & ANNOTATION ---------------------------------------------------
DefaultAssay(obj) <- "sketch"
obj[["sketch"]] <- JoinLayers(obj[["sketch"]])
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5,
                          max.cells.per.ident = 500)          # install 'presto' for speed
write.csv(markers, file.path(OUT_DIR, "clusters/cluster_markers.csv"), row.names = FALSE)

# Example marker sets (mouse skin), chosen from genes ON the Xenium Mouse 5K
# panel. Many textbook markers are NOT on that panel (e.g. Krt14, Krt5, Krt10,
# Lor, Ivl, Acta2, Lyz2, C1qa, S100a8/9, Cpa3, Dcn), so always check coverage.
markers_all <- list(
  Basal_KC          = c("Trp63", "Itgb4", "Itga6", "Lama3", "Lamb3", "Cdh3"),
  Differentiated_KC = c("Dsg1a", "Dsg1b", "Dsc1", "Evpl", "Ppl", "Zfp750", "Krt78", "Alox12b"),
  HairFollicle      = c("Sox9", "Lgr5", "Lhx2", "Krt79", "Krt71", "Krt28", "Lgr6"),
  Sebocyte          = c("Pparg", "Plin2"),
  Fibroblast        = c("Col1a1", "Col1a2", "Pdgfra", "Lum", "Sfrp2"),
  Dermal_papilla    = c("Crabp1", "Sox2", "Alx4"),
  Endothelial       = c("Pecam1", "Cdh5", "Kdr", "Plvap", "Flt1"),
  Lymphatic_EC      = c("Prox1", "Lyve1", "Flt4"),
  Pericyte_SMC      = c("Myh11", "Pdgfrb", "Notch3", "Tagln", "Kcnj8", "Abcc9"),
  Skeletal_muscle   = c("Myh4", "Myh1", "Myh2", "Des"),
  Macrophage        = c("Adgre1", "Cd68", "Mrc1", "Csf1r", "Cd163", "Ms4a7", "Fcgr1"),
  Langerhans        = c("Cd207", "Epcam"),
  Dendritic_cell    = c("Flt3", "Xcr1", "Clec9a", "Cd209a", "Itgax", "Ccr7"),
  T_cell            = c("Cd3e", "Cd3g", "Cd3d", "Cd2", "Lck", "Cd8a"),
  NK_cell           = c("Klrb1c", "Ncr1", "Gzma", "Klrk1"),
  B_cell            = c("Cd79a", "Cd79b", "Ms4a1", "Cd19"),
  Plasma_cell       = c("Jchain", "Mzb1"),
  Neutrophil        = c("Ly6g", "Csf3r", "Cxcr2"),
  Mast_cell         = c("Kit", "Tpsb2", "Ms4a2"),
  Melanocyte        = c("Dct", "Tyrp1", "Tyr", "Mitf"),
  Schwann_nerve     = c("Sox10", "S100b", "Ngfr"),
  Adipocyte         = c("Adipoq", "Plin1", "Lep", "Cidea"),
  Cycling           = c("Mki67", "Top2a", "Ccnb1"))   # a state, not a type: check which lineage cycles
panel <- rownames(obj[["sketch"]])
cov <- data.frame(set = names(markers_all), n_markers = lengths(markers_all),
                  n_on_panel = sapply(markers_all, function(g) sum(g %in% panel)),
                  missing = sapply(markers_all, function(g) paste(setdiff(g, panel), collapse = ",")))
cov$used <- cov$n_on_panel >= 2
write.csv(cov, file.path(OUT_DIR, "clusters/marker_set_coverage.csv"), row.names = FALSE)
use_markers <- lapply(markers_all[cov$used], intersect, panel)

# Score the sketch cells in a small stand-alone object (AddModuleScore fails
# on a sketch assay that holds fewer cells than the full object).
sk <- CreateSeuratObject(LayerData(obj, "counts", assay = "sketch"), assay = "sketch")
LayerData(sk, "data") <- LayerData(obj, "data", assay = "sketch")
sk <- AddModuleScore(sk, features = use_markers, name = "ct_", ctrl = min(100, floor(length(panel) / 30)))
score_cols <- paste0("ct_", seq_along(use_markers))
sc <- sk@meta.data[, score_cols]; sc$seurat_clusters <- obj$seurat_clusters[rownames(sc)]
rm(sk); gc()
cl_score <- aggregate(sc[, score_cols], list(cluster = sc$seurat_clusters), mean)
colnames(cl_score)[-1] <- names(use_markers)
sm <- as.matrix(cl_score[, names(use_markers)])
best <- apply(sm, 1, max); second <- apply(sm, 1, function(x) sort(x, decreasing = TRUE)[2])
cl_score$best_score <- round(best, 3); cl_score$margin <- round(best - second, 3)
cl_score$auto_label <- names(use_markers)[apply(sm, 1, which.max)]
cl_score$auto_label[best <= ANNOT_MIN_SCORE] <- "Unassigned"
cl_score$auto_label[best > ANNOT_MIN_SCORE & (best - second) < ANNOT_MIN_MARGIN] <- "Ambiguous"
cl_score$runner_up <- names(use_markers)[apply(sm, 1, function(x) order(x, decreasing = TRUE)[2])]
top_genes <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 8) %>%
  summarise(top_markers = paste(gene, collapse = ","), .groups = "drop")
cl_score$top_markers <- top_genes$top_markers[match(cl_score$cluster, top_genes$cluster)]
write.csv(cl_score, file.path(OUT_DIR, "clusters/cluster_module_scores.csv"), row.names = FALSE)
ggsave(file.path(OUT_DIR, "clusters/marker_dotplot.pdf"), width = 26, height = 9,
       DotPlot(obj, features = use_markers, group.by = "seurat_clusters") + RotatedAxis() +
         theme(axis.text.x = element_text(size = 6)))

# Reviewed labels: edit 'celltype' in curated/cluster_labels.csv (rename, merge
# by giving the same name, keep "Unassigned"/"Ambiguous" if unsure), set
# reviewed = TRUE and re-run. A fingerprint guards against old labels landing
# on new clusters that merely share a number.
clustering_id <- clustering_fingerprint(colnames(obj[["sketch"]]), obj$seurat_clusters[colnames(obj[["sketch"]])])
labels_file <- file.path(CURATED_DIR, "cluster_labels.csv")
if (file.exists(labels_file)) {
  old <- read.csv(labels_file, stringsAsFactors = FALSE)
  if (!"clustering_id" %in% colnames(old) || any(old$clustering_id != clustering_id)) {
    bak <- sub("\\.csv$", paste0("_OLD_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"), labels_file)
    file.rename(labels_file, bak)
    message("WARNING: clustering changed since cluster_labels.csv was written; old file saved as ",
            basename(bak), ". A new UNREVIEWED template is written.")
  }
}
if (!file.exists(labels_file))
  write.csv(data.frame(cluster = cl_score$cluster, auto_label = cl_score$auto_label, runner_up = cl_score$runner_up,
                       top_markers = cl_score$top_markers, celltype = cl_score$auto_label,
                       reviewed = FALSE, clustering_id = clustering_id), labels_file, row.names = FALSE)
labs <- read.csv(labels_file, stringsAsFactors = FALSE)
if (!setequal(as.character(labs$cluster), as.character(cl_score$cluster)))
  stop("cluster_labels.csv lists different clusters than the current clustering.")
if (any(is.na(labs$celltype) | !nzchar(trimws(labs$celltype))))
  stop("Blank celltype in cluster_labels.csv - use 'Unassigned' if unsure.")
LABELS_REVIEWED <- nrow(labs) > 0 && all(labs$reviewed %in% TRUE)
message("Cluster labels reviewed: ", LABELS_REVIEWED)
file.copy(c(core_map_file, labels_file), file.path(OUT_DIR, c("cores", "clusters")), overwrite = TRUE)
obj$celltype <- labs$celltype[match(as.character(obj$cluster_full), as.character(labs$cluster))]
obj$celltype[is.na(obj$celltype)] <- "Unassigned"
DefaultAssay(obj) <- "Xenium"

# 8a. (optional) ENRICHR CELL-TYPE LOOKUP: a hint, not evidence ----------
if (RUN_ENRICHR) {
  library(enrichR); setEnrichrSite("Enrichr")
  dbs <- grep("^PanglaoDB|^CellMarker|^Tabula_Muris", listEnrichrDbs()$libraryName, value = TRUE)
  enr <- list()
  for (cl in sort(unique(as.character(markers$cluster)))) {
    genes <- markers %>% filter(cluster == cl, p_val_adj < 0.05) %>% slice_max(avg_log2FC, n = 50) %>% pull(gene)
    if (length(genes) < 5) next
    res <- tryCatch(enrichr(toupper(genes), dbs), error = function(e) NULL)
    for (db in names(res)) if (nrow(res[[db]])) {
      top <- head(res[[db]][order(res[[db]]$Adjusted.P.value), ], 3)
      enr[[paste(cl, db)]] <- data.frame(cluster = cl, library = db, term = top$Term,
                                         adj_p = signif(top$Adjusted.P.value, 3), overlap = top$Overlap)
    }
    Sys.sleep(1)
  }
  if (length(enr)) write.csv(do.call(rbind, enr), file.path(OUT_DIR, "clusters/enrichr_celltype_hits.csv"), row.names = FALSE)
}

# 8b. (optional) RCTD REFERENCE LABELS: a second opinion -----------------
if (RUN_RCTD) {
  library(spacexr)
  ref <- readRDS(path.expand(REF_RDS))
  ref_ct <- factor(gsub("[/ ]", "_", ref@meta.data[[REF_CT_COL]])); names(ref_ct) <- colnames(ref)
  ref_cells <- names(ref_ct)[ref_ct %in% names(which(table(ref_ct) >= 25))]   # RCTD needs >= 25 cells/type
  ref_counts <- LayerData(ref, "counts", assay = "RNA")[, ref_cells]
  reference <- Reference(ref_counts, droplevels(ref_ct[ref_cells]), colSums(ref_counts)); rm(ref); gc()
  q_cells  <- if (RCTD_ALL_CELLS) colnames(obj) else colnames(obj[["sketch"]])
  q_counts <- LayerData(obj, "counts", assay = "Xenium")[, q_cells]
  q_xy <- data.frame(x = obj@meta.data[q_cells, "x_centroid"], y = obj@meta.data[q_cells, "y_centroid"], row.names = q_cells)
  rctd <- create.RCTD(SpatialRNA(q_xy, q_counts, colSums(q_counts)), reference, max_cores = 8,
                      UMI_min = 20, counts_MIN = 10, UMI_min_sigma = 100)     # Xenium cells are low-count
  rctd <- run.RCTD(rctd, doublet_mode = "doublet")
  res  <- rctd@results$results_df
  ft <- as.character(res$first_type); ft[res$spot_class == "reject"] <- NA   # rejected = no call
  obj$rctd_celltype <- NA_character_; obj$rctd_celltype[match(rownames(res), colnames(obj))] <- ft
  rc <- obj@meta.data[!is.na(obj$rctd_celltype), ]
  write.csv(as.data.frame.matrix(table(cluster = rc$cluster_full, rctd = rc$rctd_celltype)),
            file.path(OUT_DIR, "clusters/rctd_by_cluster.csv"))
  saveRDS(rctd, file.path(OUT_DIR, "clusters/rctd_result.rds"))
}

# 9. SPATIAL PLOTS --------------------------------------------------------
ct_lv <- sort(unique(obj$celltype))
ct_cols <- setNames(DiscretePalette(max(3, length(ct_lv)), "polychrome")[seq_along(ct_lv)], ct_lv)
genes_show <- intersect(GENES_SHOW, rownames(obj))
for (s in names(dirs)) {
  ggsave(file.path(OUT_DIR, "spatial", paste0(s, "_celltypes.pdf")), width = 12, height = 12,
         ImageDimPlot(obj, fov = s, group.by = "celltype", cols = ct_cols, size = 0.3, dark.background = FALSE))
  if (length(genes_show))
    ggsave(file.path(OUT_DIR, "spatial", paste0(s, "_genes.pdf")), width = 18, height = 12,
           ImageFeaturePlot(obj, fov = s, features = genes_show, size = 0.3, max.cutoff = "q95"))
}

# 9b. (optional) ZOOMED CORES WITH CELL OUTLINES --------------------------
# Outlines are read per core straight from cell_boundaries.parquet (arrow
# filters on disk), so the whole slide's polygons are never loaded.
if (RUN_ZOOM) {
  pdf(file.path(OUT_DIR, "spatial/core_zooms.pdf"), width = 12, height = 10)
  for (s in names(dirs)) {
    d <- obj@meta.data[obj$slide == s, ]; d$cell_id <- sub(paste0("^", s, "_"), "", rownames(d))
    bnd_ds <- arrow::open_dataset(file.path(dirs[[s]], "cell_boundaries.parquet"))
    for (cid in head(names(sort(table(d$core_id), decreasing = TRUE)), ZOOM_N)) {
      dc <- d[d$core_id == cid, ]; ids <- dc$cell_id
      bnd <- bnd_ds |> dplyr::filter(cell_id %in% ids) |> dplyr::collect() |> as.data.frame()
      bnd$cell_id <- as.character(bnd$cell_id); bnd$celltype <- dc$celltype[match(bnd$cell_id, dc$cell_id)]
      print(ggplot(bnd, aes(vertex_x, vertex_y, group = cell_id, fill = celltype)) +
              geom_polygon(colour = "white", linewidth = 0.05) + scale_fill_manual(values = ct_cols) +
              scale_y_reverse() + coord_fixed() + theme_void() +
              theme(panel.background = element_rect(fill = "black")) +
              ggtitle(paste(s, cid, dc$mouse_id[1], dc$condition[1], sep = " | ")))
    }
  }
  dev.off()
}

# 10. CELLULAR NICHES (within cores, clustered jointly over all cores) ----
md <- obj@meta.data
np <- niche_profiles(md, sort(unique(md$celltype)), NICHE_K, NICHE_MAX_DIST)
obj$n_neighbours <- np$n_neighbours[colnames(obj)]
message("Isolated cells (no neighbour within ", NICHE_MAX_DIST, " um): ", sum(obj$n_neighbours == 0))
if (is.null(np$mat) || nrow(unique(round(np$mat, 3))) < 3) {
  message("WARNING: too few distinct neighbourhoods - niche analysis skipped."); obj$niche <- factor("Isolated")
} else {
  k_use <- min(N_NICHES, nrow(unique(round(np$mat, 3))) - 1)
  km <- kmeans(np$mat, centers = k_use, nstart = 5, iter.max = 50)
  niche_vec <- setNames(paste0("N", km$cluster), rownames(np$mat))
  obj$niche <- unname(niche_vec[colnames(obj)]); obj$niche[is.na(obj$niche)] <- "Isolated"
  obj$niche <- factor(obj$niche)
  niche_comp <- aggregate(np$mat, list(niche = niche_vec[rownames(np$mat)]), mean)
  write.csv(niche_comp, file.path(OUT_DIR, "niches/niche_celltype_composition.csv"), row.names = FALSE)
  ggsave(file.path(OUT_DIR, "niches/niche_heatmap.pdf"), width = 11, height = 5,
         ggplot(pivot_longer(niche_comp, -niche, names_to = "celltype", values_to = "fraction"),
                aes(celltype, niche, fill = fraction)) + geom_tile() + scale_fill_viridis_c() + theme_bw() + RotatedAxis())
}

# 10b. (optional) BANKSY SPATIAL DOMAINS ---------------------------------
# Cores are shifted 1e6 um apart before the neighbour search, so Banksy
# neighbours never come from another core; very small cores are skipped.
if (RUN_BANKSY) {
  library(Banksy); library(SpatialExperiment); library(scuttle)
  hvg <- VariableFeatures(obj[["sketch"]]); obj$banksy_domain <- NA_character_
  for (s in names(dirs)) {
    core_n <- table(obj$sample_id[obj$slide == s]); small <- names(core_n)[core_n < 31]
    cells <- colnames(obj)[obj$slide == s & !obj$sample_id %in% small]
    xy <- as.matrix(obj@meta.data[cells, c("x_centroid", "y_centroid")])
    xy[, 1] <- xy[, 1] + 1e6 * as.integer(factor(obj@meta.data[cells, "sample_id"]))
    se <- SpatialExperiment(assays = list(counts = LayerData(obj, "counts", assay = "Xenium")[hvg, cells]),
                            spatialCoords = xy)
    se <- computeLibraryFactors(se); assay(se, "normcounts") <- normalizeCounts(se, log = FALSE)
    se <- computeBanksy(se, assay_name = "normcounts", compute_agf = TRUE, k_geom = 15)
    pca_args <- list(se, use_agf = TRUE, lambda = 0.8, npcs = 20)
    if ("lazy" %in% names(formals(runBanksyPCA))) pca_args$lazy <- FALSE   # newer Banksy
    se <- do.call(runBanksyPCA, pca_args)
    se <- clusterBanksy(se, use_agf = TRUE, lambda = 0.8, resolution = 0.8, npcs = 20)
    obj$banksy_domain[match(cells, colnames(obj))] <- paste0(s, "_D", colData(se)[[clusterNames(se)[1]]])
    rm(se); gc()
  }
}

# 10c. EXPORT TO XENIUM EXPLORER & SAVE ------------------------------------
for (s in names(dirs)) {
  d <- obj@meta.data[obj$slide == s, ]; ids <- sub(paste0("^", s, "_"), "", rownames(d))
  for (g in intersect(c("cluster_full", "celltype", "niche", "core_id", "rctd_celltype", "banksy_domain"), colnames(d)))
    write.csv(data.frame(cell_id = ids, group = d[[g]]),
              file.path(OUT_DIR, "explorer", paste0(s, "_", g, ".csv")), row.names = FALSE, quote = FALSE)
}
saveRDS(obj, file.path(OUT_DIR, "xenium_tma_seurat.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

# 10d. SAMPLE TREE: animal > condition > core > cells -----------------------
tree <- obj@meta.data %>% group_by(condition, stratum, mouse_id, sample_id, slide) %>%
  summarise(cells = n(), median_counts = median(nCount_Xenium), .groups = "drop")
write.csv(tree, file.path(OUT_DIR, "sample_tree.csv"), row.names = FALSE)
ggsave(file.path(OUT_DIR, "sample_tree.pdf"), width = 12, height = 5,
       ggplot(tree, aes(mouse_id, cells, fill = slide, group = sample_id)) +
         geom_col(position = position_dodge2(preserve = "single"), colour = "white") +
         geom_hline(yintercept = 1000, linetype = "dashed", colour = "grey50") +
         facet_grid(~ condition, scales = "free_x", space = "free_x") +
         labs(x = NULL, y = "Cells per core (after QC)",
              subtitle = "bars side by side = cores of one animal (not separate replicates)") +
         theme_classic() + theme(axis.text.x = element_text(angle = 60, hjust = 1)))

# 11. COMPOSITION TABLES (always; per animal when metadata are complete) ---
md <- obj@meta.data
md$rep_id <- if (META_OK) paste(md$mouse_id, md$condition, sep = "|") else md$sample_id
prop_ct <- md %>% count(rep_id, condition, celltype) %>% group_by(rep_id) %>% mutate(frac = n / sum(n)) %>% ungroup()
tag <- if (META_OK) "per_animal" else "per_core_EXPLORATORY"
write.csv(prop_ct, file.path(OUT_DIR, "de", paste0("celltype_proportions_", tag, ".csv")), row.names = FALSE)

# 12. STATISTICS: composition test + pseudobulk DE ------------------------
if (!(META_OK && LABELS_REVIEWED)) {
  message("STATISTICS SKIPPED: need complete core_map.csv (META_OK = ", META_OK,
          ") and reviewed cluster_labels.csv (LABELS_REVIEWED = ", LABELS_REVIEWED, ").")
} else {
  for (p in c("DESeq2", "limma")) if (!requireNamespace(p, quietly = TRUE)) stop("Install ", p)
  lv <- c(COMPARE_REF, COMPARE_TEST)
  md2 <- md
  if (!is.na(DE_SUBSET_STRATUM)) md2 <- md2[md2$stratum == DE_SUBSET_STRATUM, ]
  if (!all(lv %in% md2$condition)) stop("Groups ", paste(lv, collapse = " / "), " not found in condition")
  md2 <- md2[md2$condition %in% lv, ]; md2$grp <- factor(md2$condition, levels = lv)
  if (!PAIRED && any(tapply(md2$grp, md2$mouse_id, function(x) length(unique(x))) > 1))
    stop("Some animals appear in both groups but PAIRED = FALSE.")
  gs <- table(md2$grp, md2$slide) > 0
  if (all(rowSums(gs) == 1) && all(colSums(gs) <= 1))
    message("NOTE: each group sits on its own slide - group and slide effects are inseparable; ",
            "treat the results below as exploratory.")
  cname <- paste0(lv[2], "_vs_", lv[1], "__stratum_", ifelse(is.na(DE_SUBSET_STRATUM), "all", DE_SUBSET_STRATUM),
                  "__", ifelse(PAIRED, "paired", "unpaired"))
  STATS_DIR <- file.path(OUT_DIR, "de", cname); dir.create(STATS_DIR, recursive = TRUE, showWarnings = FALSE)
  comp <- composition_test(md2, lv, PAIRED, MIN_MICE)
  if (!is.null(comp)) write.csv(comp, file.path(STATS_DIR, "composition_test.csv"), row.names = FALSE)
  pbde <- pseudobulk_de(LayerData(obj, "counts", assay = "Xenium"), md2, lv, PAIRED, MIN_MICE, MIN_CELLS_PB)
  for (ct in names(pbde$results))
    write.csv(pbde$results[[ct]], file.path(STATS_DIR, paste0("pseudobulk_DE_", gsub("[^A-Za-z0-9]", "_", ct), ".csv")),
              row.names = FALSE)
  write.csv(pbde$info, file.path(STATS_DIR, "pseudobulk_samples.csv"), row.names = FALSE)
  if (!is.null(pbde$summary)) write.csv(pbde$summary, file.path(STATS_DIR, "DE_summary.csv"), row.names = FALSE)
  if (!is.null(pbde$skipped)) write.csv(pbde$skipped, file.path(STATS_DIR, "DE_skipped.csv"), row.names = FALSE)
  message("DE (", cname, "): ", length(pbde$results), " cell types tested; see DE_skipped.csv for the rest.")
}
message("Done. Results in ", OUT_DIR)

# ======================================================================
# NEXT STEPS (in order)
# ======================================================================
# Outputs of a run: runs/<RUN_NAME>/ (path also in latest_run.txt).
# Hand-edited inputs: curated/ (read by every run).
# 1. Cores: open cores/*_core_map.pdf and *_core_map_matched.pdf and compare
#    a few cores with the H&E / TMA map. *_orientations.csv lists all 8
#    orientations; the chosen one should have 0 mismatches. Split core ->
#    raise CORE_EPS / CORE_MERGE_DIST; merged cores -> lower them; lost small
#    core -> lower CORE_MIN_CELLS. After changes, delete curated/core_map.csv.
# 2. Metadata: check curated/core_map.csv; every kept core needs mouse_id,
#    condition and stratum; set exclude = TRUE for folded / damaged cores.
#    Look at cores/condition_by_slide.csv: every group should be on both
#    slides. If not, group and slide effects cannot be separated.
# 3. QC: qc/qc_per_core.csv says why each core lost cells. Fix causes (folds,
#    segmentation) before changing thresholds; keep one rule for all samples.
# 4. Batch check: clusters/umap.pdf by slide; the same cell types should overlap.
# 5. Labels: compare cluster_module_scores.csv, marker_dotplot.pdf,
#    cluster_markers.csv, spatial plots (and RCTD / Enrichr if run). Edit
#    curated/cluster_labels.csv, keep doubtful clusters Unassigned/Ambiguous,
#    set reviewed = TRUE, re-run (same settings -> same clustering, labels kept).
# 6. Explorer: load explorer/*.csv (Cell > Add cell categorization) and check
#    labels on DAPI / H&E.
# 7. Niches: name niches from niches/niche_heatmap.pdf; try NICHE_K 15 / 50 and
#    N_NICHES 6-12; only name niches that stay stable.
# 8. Statistics run once steps 2 and 5 are done: de/<test>_vs_<ref>__.../
#    composition_test.csv, DE_summary.csv, DE_skipped.csv. >= 3 animals (or
#    complete pairs) per group recommended; report effect sizes.
# 9. Pathways: use the PANEL genes (rownames(obj)) as background, not the genome.
# 10. Methods: each run folder keeps run_settings.txt, sessionInfo.txt and the
#     core_map.csv / cluster_labels.csv it used.
