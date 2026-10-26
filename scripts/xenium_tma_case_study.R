# ======================================================================
# XENIUM TMA CASE STUDY – multi-slide tissue microarray pipeline (R)
# OmicLearningHub · Module 1 · https://githubrudramani.github.io/OmicLearningHub/xenium.html
#
# Example design: two Xenium slides, each a tissue microarray (TMA) of
# mouse tissue cores; each slide holds one model (model_1, model_2), and
# each mouse contributes one or more cores. Marker lists below use mouse
# skin as the example tissue – replace them with markers for your tissue.
#
# What this script shows beyond the basic tutorial:
#   * automatic TMA core detection (DBSCAN) with a core map PDF to check it
#   * mouse (not core, not cell) as the statistical replicate
#   * Seurat sketching for >1M cells
#   * no Harmony when slide = biological group (it would remove the biology)
#   * marker-module annotation, optional RCTD / Banksy, niches,
#     per-mouse composition and pseudobulk DE
# ======================================================================

library(Seurat); library(arrow); library(dbscan)
library(ggplot2); library(dplyr); library(Matrix)
options(future.globals.maxSize = 50 * 1024^3)
set.seed(1)

# Optional steps (FALSE = skip; switch on one at a time after the first run)
RUN_RCTD   <- FALSE   # 8b: reference-based labels (package spacexr)
REF_RDS    <- "path/to/your_scRNAseq_reference.rds"   # build one with scripts/build_skin_rctd_reference.R
REF_CT_COL <- "celltype"   # the cell-type column in that object
RCTD_ALL_CELLS <- FALSE  # FALSE = sketch cells only (fast); TRUE = every cell (hours)
RUN_ZOOM   <- TRUE    # 9b: zoomed cores with cell outlines (Seurat Crop)
ZOOM_N     <- 4       # largest cores per slide to zoom
RUN_BANKSY <- FALSE   # 10b: Banksy spatial domains (memory-heavy)

# Paths: point these at your two Xenium output bundles
DIR_SLIDE1 <- path.expand("path/to/xenium_output_slide1")
DIR_SLIDE2 <- path.expand("path/to/xenium_output_slide2")
OUT_DIR    <- path.expand("xenium_tma_results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (d in c("qc", "cores", "clusters", "spatial", "niches", "de", "explorer"))
  dir.create(file.path(OUT_DIR, d), showWarnings = FALSE)

# 1. LOAD SLIDES ---------------------------------------------------------
load_slide <- function(dir, slide, model) {
  o <- LoadXenium(data.dir = dir, fov = slide, assay = "Xenium", molecule.coordinates = FALSE)
  o$slide <- slide; o$model <- model
  cm <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  rownames(cm) <- as.character(cm$cell_id)
  AddMetaData(o, cm[colnames(o), c("x_centroid", "y_centroid", "cell_area",
                                   "control_probe_counts", "control_codeword_counts")])
}
tma1 <- load_slide(DIR_SLIDE1, "TMA1", "model_1")
tma2 <- load_slide(DIR_SLIDE2, "TMA2", "model_2")

# 2. MERGE ---------------------------------------------------------------
obj <- merge(tma1, y = tma2, add.cell.ids = c("TMA1", "TMA2"))
rm(tma1, tma2); gc()
DefaultAssay(obj) <- "Xenium"
obj[["Xenium"]] <- JoinLayers(obj[["Xenium"]])
obj$neg_frac <- rowSums(cbind(obj$control_probe_counts, obj$control_codeword_counts),
                        na.rm = TRUE) / pmax(obj$nCount_Xenium, 1)

# 3. DETECT TMA CORES (DBSCAN -> merge fragments -> drop debris) ---------
# Note: match cells with colnames(obj); rownames(obj) are genes.
obj$core_id <- NA_character_
for (s in c("TMA1", "TMA2")) {
  cells <- colnames(obj)[obj$slide == s]
  xy <- as.matrix(obj@meta.data[cells, c("x_centroid", "y_centroid")])
  cl <- dbscan::dbscan(xy, eps = 75, minPts = 10)$cluster
  cen <- aggregate(xy[cl > 0, ], list(cl = cl[cl > 0]), median)
  grp <- cutree(hclust(dist(cen[, -1]), method = "single"), h = 400)  # merge pieces < 400 um apart
  cl[cl > 0] <- grp[match(cl[cl > 0], cen$cl)]
  big <- as.integer(names(which(table(cl[cl > 0]) >= 300)))           # drop clusters < 300 cells
  obj$core_id[match(cells, colnames(obj))] <- ifelse(cl %in% big, paste0("C", cl), NA)

  d <- obj@meta.data[cells, ]
  lab <- d %>% filter(!is.na(core_id)) %>% group_by(core_id) %>%
    summarise(x = median(x_centroid), y = median(y_centroid), n = n())
  p <- ggplot(d[sample(nrow(d), min(nrow(d), 2e5)), ],
              aes(x_centroid, y_centroid, colour = is.na(core_id))) +
    geom_point(size = 0.05) + scale_y_reverse() + coord_fixed() +
    scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "red"), guide = "none") +
    geom_label(data = lab, aes(x, y, label = core_id), inherit.aes = FALSE, size = 2.5) +
    theme_void() + ggtitle(paste(s, "- detected cores (red = unassigned)"))
  ggsave(file.path(OUT_DIR, "cores", paste0(s, "_core_map.pdf")), p, width = 12, height = 12)
}
obj$sample_id <- paste(obj$slide, obj$core_id, sep = "_")

# Core -> mouse mapping. The replicate for statistics is the MOUSE: several
# cores from one mouse are not independent. Fill cores/core_map.csv
# (columns: sample_id, mouse_id, model, exclude) from your TMA layout sheet
# and the core map PDFs; until then each core counts as its own mouse.
core_map_file <- file.path(OUT_DIR, "cores", "core_map.csv")
if (!file.exists(core_map_file)) {
  tmpl <- obj@meta.data %>% filter(!is.na(core_id)) %>% count(slide, core_id, sample_id, model) %>%
    mutate(mouse_id = sample_id, exclude = FALSE)
  write.csv(tmpl, core_map_file, row.names = FALSE)
}
core_map <- read.csv(core_map_file, stringsAsFactors = FALSE)
k <- match(obj$sample_id, core_map$sample_id)
obj$mouse_id <- core_map$mouse_id[k]
obj$model    <- ifelse(is.na(k), obj$model, core_map$model[k])
obj$exclude  <- !is.na(k) & as.logical(core_map$exclude[k])

# 4. QC --------------------------------------------------------------------
max_area <- quantile(obj$cell_area, 0.995, na.rm = TRUE)
keep <- colnames(obj)[obj$nCount_Xenium >= 20 & obj$nFeature_Xenium >= 15 &
                      obj$neg_frac <= 0.05 & obj$cell_area >= 10 &
                      obj$cell_area <= max_area & !is.na(obj$core_id) & !obj$exclude]
obj <- subset(obj, cells = keep)
print(table(obj$slide))
obj[["Xenium"]] <- split(obj[["Xenium"]], f = obj$slide)

# 5. NORMALIZE & SKETCH (25k cells per slide) ---------------------------
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- SketchData(obj, ncells = 25000, method = "LeverageScore", sketched.assay = "sketch")

# 6. CLUSTER THE SKETCH (no Harmony: slide = model, it would remove biology)
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, return.model = TRUE, reduction.name = "umap")

# 7. PROJECT TO ALL CELLS -------------------------------------------------
obj <- ProjectData(obj, assay = "Xenium", sketched.assay = "sketch",
                   sketched.reduction = "pca", full.reduction = "pca.full",
                   dims = 1:30, umap.model = "umap",
                   refdata = list(cluster_full = "seurat_clusters"))
DefaultAssay(obj) <- "Xenium"
obj[["Xenium"]] <- JoinLayers(obj[["Xenium"]])
obj$cluster_full <- factor(obj$cluster_full, levels = sort(as.numeric(unique(obj$cluster_full))))
ggsave(file.path(OUT_DIR, "clusters/umap.pdf"),
       DimPlot(obj, reduction = "full.umap", group.by = c("cluster_full", "slide"),
               label = TRUE, raster = TRUE), width = 14, height = 6)

# 8. MARKERS & AUTOMATED ANNOTATION ---------------------------------------
DefaultAssay(obj) <- "sketch"
obj[["sketch"]] <- JoinLayers(obj[["sketch"]])
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5,
                          max.cells.per.ident = 500)
write.csv(markers, file.path(OUT_DIR, "clusters/cluster_markers.csv"), row.names = FALSE)

# Example marker sets (mouse skin). Replace with markers for your tissue.
tissue_markers <- list(
  Basal_KC        = c("Krt14", "Krt5", "Krt15", "Col17a1", "Trp63"),
  Suprabasal_KC   = c("Krt10", "Krt1", "Lor", "Ivl", "Flg"),
  HairFollicle    = c("Krt17", "Sox9", "Lgr5", "Krt79", "Krt25", "Krt71", "Lhx2"),
  Sebocyte        = c("Scd1", "Mgst1", "Pparg", "Elovl3", "Awat1"),
  Fibroblast      = c("Col1a1", "Col1a2", "Pdgfra", "Dcn", "Lum"),
  Dermal_papilla  = c("Crabp1", "Corin", "Sox2", "Alx4"),
  Endothelial     = c("Pecam1", "Cdh5", "Kdr", "Emcn"),
  Lymphatic_EC    = c("Prox1", "Lyve1", "Flt4"),
  Pericyte_SMC    = c("Acta2", "Rgs5", "Myh11", "Pdgfrb", "Des"),
  Macrophage      = c("Adgre1", "Cd68", "Lyz2", "C1qa", "Mrc1"),
  Langerhans_DC   = c("Cd207", "Epcam", "H2-Ab1", "Itgax"),
  T_cell          = c("Cd3e", "Cd3g", "Trdc", "Cd8a", "Cd4"),
  Neutrophil      = c("S100a8", "S100a9", "Ly6g", "Csf3r"),
  Mast_cell       = c("Cpa3", "Mcpt4", "Kit", "Cma1"),
  Melanocyte      = c("Pmel", "Dct", "Tyrp1", "Mlana"),
  Schwann_nerve   = c("Mpz", "Plp1", "Sox10", "S100b"),
  Adipocyte       = c("Adipoq", "Plin1", "Lep"),
  Skeletal_muscle = c("Acta1", "Myh4", "Ckm"))
tissue_markers <- lapply(tissue_markers, intersect, rownames(obj[["sketch"]]))
tissue_markers <- tissue_markers[lengths(tissue_markers) >= 2]   # keep sets with >= 2 panel genes

obj <- AddModuleScore(obj, features = tissue_markers, name = "ct_", assay = "sketch")
score_cols <- paste0("ct_", seq_along(tissue_markers))
sc <- obj@meta.data[!is.na(obj$seurat_clusters), c("seurat_clusters", score_cols)]
cl_score <- aggregate(sc[, score_cols], list(cluster = sc$seurat_clusters), mean)
colnames(cl_score)[-1] <- names(tissue_markers)
cl_score$auto_label <- names(tissue_markers)[apply(cl_score[, -1], 1, which.max)]
write.csv(cl_score, file.path(OUT_DIR, "clusters/cluster_module_scores.csv"), row.names = FALSE)
# Check these auto labels against cluster_markers.csv before trusting them.
obj$celltype <- cl_score$auto_label[match(as.character(obj$cluster_full), as.character(cl_score$cluster))]
DefaultAssay(obj) <- "Xenium"

# 8b. (optional) RCTD REFERENCE LABELS ------------------------------------
# Needs a matched scRNA-seq Seurat object (same species and tissue).
# Result: obj$rctd_celltype, plus a table comparing it with obj$celltype.
if (RUN_RCTD) {
  library(spacexr)
  ref <- readRDS(path.expand(REF_RDS))
  ref_ct <- factor(gsub("[/ ]", "_", ref@meta.data[[REF_CT_COL]]))   # RCTD dislikes "/"
  names(ref_ct) <- colnames(ref)
  keep_ct <- names(which(table(ref_ct) >= 25))                         # RCTD needs >= 25 cells/type
  ref_cells <- names(ref_ct)[ref_ct %in% keep_ct]
  ref_counts <- LayerData(ref, "counts", assay = "RNA")[, ref_cells]
  reference <- Reference(ref_counts, droplevels(ref_ct[ref_cells]), colSums(ref_counts))
  rm(ref); gc()

  q_cells <- if (RCTD_ALL_CELLS) colnames(obj) else colnames(obj[["sketch"]])
  q_counts <- LayerData(obj, "counts", assay = "Xenium")[, q_cells]
  q_coords <- data.frame(x = obj@meta.data[q_cells, "x_centroid"],
                         y = obj@meta.data[q_cells, "y_centroid"], row.names = q_cells)
  query <- SpatialRNA(q_coords, q_counts, colSums(q_counts))
  rctd  <- create.RCTD(query, reference, max_cores = 8,
                       UMI_min = 20, counts_MIN = 10, UMI_min_sigma = 100)  # Xenium cells are low-count
  rctd  <- run.RCTD(rctd, doublet_mode = "doublet")
  res <- rctd@results$results_df
  obj$rctd_celltype <- NA_character_
  obj$rctd_celltype[match(rownames(res), colnames(obj))] <- as.character(res$first_type)
  obj$rctd_class    <- NA_character_
  obj$rctd_class[match(rownames(res), colnames(obj))] <- as.character(res$spot_class)
  write.csv(as.data.frame.matrix(table(marker = obj$celltype, rctd = obj$rctd_celltype)),
            file.path(OUT_DIR, "clusters/marker_vs_rctd_labels.csv"))
  saveRDS(rctd, file.path(OUT_DIR, "clusters/rctd_result.rds"))
}

# 9. SPATIAL PLOTS --------------------------------------------------------
for (s in c("TMA1", "TMA2")) {
  ggsave(file.path(OUT_DIR, "spatial", paste0(s, "_celltypes.pdf")),
         ImageDimPlot(obj, fov = s, group.by = "celltype", size = 0.3,
                      dark.background = FALSE) + ggtitle(s), width = 12, height = 12)
}

# 9b. (optional) ZOOM INTO CORES WITH CELL OUTLINES (Seurat Crop) --------
# Crop box = the core's own cell centroids (microns), so no hand-typed coords.
if (RUN_ZOOM) {
  ct_cols <- setNames(DiscretePalette(length(unique(obj$celltype)), "polychrome"),
                      sort(unique(obj$celltype)))
  pdf(file.path(OUT_DIR, "spatial/core_zooms.pdf"), width = 12, height = 10)
  for (s in c("TMA1", "TMA2")) {
    d <- obj@meta.data[obj$slide == s, ]
    for (cid in head(names(sort(table(d$core_id), decreasing = TRUE)), ZOOM_N)) {
      dc <- d[d$core_id == cid, ]
      zoom <- Crop(obj[[s]], x = range(dc$x_centroid) + c(-20, 20),
                   y = range(dc$y_centroid) + c(-20, 20), coords = "tissue")
      if ("segmentation" %in% Boundaries(zoom)) DefaultBoundary(zoom) <- "segmentation"
      obj[["zoom"]] <- zoom
      print(ImageDimPlot(obj, fov = "zoom", group.by = "celltype", cols = ct_cols,
                         border.color = "white", border.size = 0.1, axes = TRUE) +
              ggtitle(paste(s, cid, dc$mouse_id[1], dc$model[1], sep = " | ")))
      obj[["zoom"]] <- NULL
    }
  }
  dev.off()
}

# 10. CELLULAR NICHES (joint k-means over both slides -> comparable) -----
md <- obj@meta.data
ct_levels <- sort(unique(md$celltype))
nb_list <- list()
for (s in c("TMA1", "TMA2")) {
  d   <- md[md$slide == s, ]
  knn <- dbscan::kNN(as.matrix(d[, c("x_centroid", "y_centroid")]), k = 30)$id
  ct  <- factor(d$celltype, levels = ct_levels)
  m   <- sparseMatrix(i = rep(seq_len(nrow(d)), 30), j = as.integer(ct[as.vector(knn)]), x = 1,
                      dims = c(nrow(d), length(ct_levels)), dimnames = list(rownames(d), ct_levels))
  nb_list[[s]] <- as.matrix(m) / 30
}
nb_mat <- do.call(rbind, nb_list)
km <- kmeans(nb_mat, centers = 8, nstart = 5, iter.max = 50)
obj$niche <- factor(paste0("N", km$cluster[match(colnames(obj), rownames(nb_mat))]))
write.csv(aggregate(nb_mat, list(niche = paste0("N", km$cluster)), mean),
          file.path(OUT_DIR, "niches/niche_celltype_composition.csv"), row.names = FALSE)

# 10b. (optional) BANKSY SPATIAL DOMAINS ----------------------------------
# Clusters cells on their own + neighbours' expression. lambda 0.8 = tissue
# domains; 0.2 = closer to cell typing.
# Runs per slide, so domain numbers are not shared between slides.
if (RUN_BANKSY) {
  library(Banksy); library(SpatialExperiment); library(scuttle)
  hvg <- VariableFeatures(obj[["sketch"]])
  obj$banksy_domain <- NA_character_
  for (s in c("TMA1", "TMA2")) {
    cells <- colnames(obj)[obj$slide == s]
    se <- SpatialExperiment(
      assays = list(counts = LayerData(obj, "counts", assay = "Xenium")[hvg, cells]),
      spatialCoords = as.matrix(obj@meta.data[cells, c("x_centroid", "y_centroid")]))
    se <- computeLibraryFactors(se)
    assay(se, "normcounts") <- normalizeCounts(se, log = FALSE)
    se <- computeBanksy(se, assay_name = "normcounts", compute_agf = TRUE, k_geom = 15)
    se <- runBanksyPCA(se, use_agf = TRUE, lambda = 0.8, npcs = 20)
    se <- clusterBanksy(se, use_agf = TRUE, lambda = 0.8, resolution = 0.8, npcs = 20)
    obj$banksy_domain[match(cells, colnames(obj))] <-
      paste0(s, "_D", colData(se)[[clusterNames(se)[1]]])
    rm(se); gc()
    ggsave(file.path(OUT_DIR, "niches", paste0(s, "_banksy_domains.pdf")),
           ImageDimPlot(obj, fov = s, group.by = "banksy_domain", size = 0.3,
                        dark.background = FALSE), width = 12, height = 12)
  }
}

# 11. COMPOSITION PER MOUSE -------------------------------------------------
md <- obj@meta.data
prop <- md %>% count(model, mouse_id, celltype) %>%
  group_by(mouse_id) %>% mutate(frac = n / sum(n)) %>% ungroup()
write.csv(prop, file.path(OUT_DIR, "de/celltype_proportions_per_mouse.csv"), row.names = FALSE)

# 12. PSEUDOBULK DE (model_2 vs model_1, replicate = mouse) ---------------
n_per_model <- tapply(md$mouse_id, md$model, function(x) length(unique(x)))
print(n_per_model)
if (all(n_per_model >= 2) && requireNamespace("DESeq2", quietly = TRUE)) {
  counts <- LayerData(obj, "counts", assay = "Xenium")
  grp <- paste(md$celltype, md$mouse_id, sep = "||")
  mm  <- sparse.model.matrix(~ 0 + factor(grp))
  colnames(mm) <- levels(factor(grp))
  pb  <- as.matrix(counts %*% mm)
  info <- data.frame(id = colnames(pb),
                     celltype = sub("\\|\\|.*", "", colnames(pb)),
                     mouse_id = sub(".*\\|\\|", "", colnames(pb)),
                     n_cells  = as.vector(table(grp)[colnames(pb)]))
  info$model <- md$model[match(info$mouse_id, md$mouse_id)]
  models <- c("model_1", "model_2")
  for (ct in unique(info$celltype)) {
    ii <- info[info$celltype == ct & info$n_cells >= 20, ]
    if (any(table(factor(ii$model, levels = models)) < 2)) next
    ii$model <- factor(ii$model, levels = models)
    dds <- DESeq2::DESeqDataSetFromMatrix(round(pb[, ii$id]), ii, ~ model)
    dds <- dds[rowSums(DESeq2::counts(dds) >= 5) >= 2, ]
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    r <- as.data.frame(DESeq2::results(dds, contrast = c("model", models[2], models[1])))
    r$gene <- rownames(r)
    write.csv(r[order(r$padj), ], file.path(OUT_DIR, "de",
              paste0("pseudobulk_DE_", gsub("[^A-Za-z0-9]", "_", ct), "_model_2_vs_model_1.csv")),
              row.names = FALSE)
  }
}

# 13. EXPORT TO XENIUM EXPLORER & SAVE -----------------------------------
for (s in c("TMA1", "TMA2")) {
  d <- obj@meta.data[obj$slide == s, ]
  ids <- sub(paste0("^", s, "_"), "", rownames(d))
  write.csv(data.frame(cell_id = ids, group = d$celltype),
            file.path(OUT_DIR, "explorer", paste0(s, "_celltype.csv")), row.names = FALSE, quote = FALSE)
  for (g in intersect(c("niche", "rctd_celltype", "banksy_domain"), colnames(d)))
    write.csv(data.frame(cell_id = ids, group = d[[g]]),
              file.path(OUT_DIR, "explorer", paste0(s, "_", g, ".csv")), row.names = FALSE, quote = FALSE)
}
saveRDS(obj, file.path(OUT_DIR, "xenium_tma_seurat.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Pipeline completed. Results saved to: ", OUT_DIR)

# ======================================================================
# NEXT STEPS (in order)
# ======================================================================
# 1. Check core detection: open cores/TMA1_core_map.pdf and TMA2_core_map.pdf.
#    Every tissue core should have exactly one label. If a core is split,
#    raise eps / the 400 um merge distance; if two cores are merged, lower it.
# 2. Fill cores/core_map.csv from your TMA layout sheet: give each core its
#    mouse_id, the real model name (e.g. WT / KO), and exclude = TRUE for
#    folded, torn or empty cores. Add a condition column if the TMA has
#    treated vs. control cores. Then re-run the script from the top.
# 3. Look at the QC numbers: print(table(obj$slide)) and the cells kept per
#    core. If one slide or core loses far more cells, loosen the thresholds
#    (counts >= 20, genes >= 15) rather than lose a whole sample.
# 4. Check for a slide batch effect: in clusters/umap.pdf, coloured by slide,
#    the same cell types should overlap. If whole clusters are made of one
#    slide only, decide whether that is biology (model) or technical before
#    adding Harmony – remember slide = model here.
# 5. Review the cell-type labels: compare clusters/cluster_module_scores.csv
#    with the top genes in clusters/cluster_markers.csv. Rename, merge or
#    split clusters by hand (or re-cluster one lineage at higher resolution)
#    before any comparison.
# 6. Validate in Xenium Explorer: load explorer/TMA*_celltype.csv and
#    TMA*_niche.csv (Cell > Add cell categorization) and check the labels sit
#    where expected on the H&E / DAPI.
# 7. Optional steps built into this script (switches at the top):
#    a. RUN_RCTD: set REF_RDS to a matched scRNA-seq Seurat object and
#       REF_CT_COL. Compare labels in clusters/marker_vs_rctd_labels.csv;
#       where they disagree, check the markers. Sketch cells first.
#    b. RUN_ZOOM: check spatial/core_zooms.pdf – are cells segmented
#       sensibly and labels in the right tissue layer?
#    c. RUN_BANKSY: if niches look noisy, compare niches/*_banksy_domains.pdf
#       with the tissue layers. Try on one slide first (memory).
# 8. Interpret niches: read niches/niche_celltype_composition.csv, name each
#    niche and compare niche fractions per mouse between the two models.
# 9. Composition test: compare cell-type proportions between models with
#    one value per mouse (e.g. speckle::propeller), not per cell.
# 10. Differential expression: read de/pseudobulk_DE_*.csv (padj < 0.05).
#     Only trust cell types with enough mice per model (>= 3 is better than 2).
#     Follow up with pathway enrichment and plot key genes spatially.
# 11. Align the H&E images to the Xenium image in Xenium Explorer
#     (Image alignment) to tie findings to histology.
# 12. Record thresholds, package versions (sessionInfo.txt) and any manual
#     label changes for the methods section.
