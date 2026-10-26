# ======================================================================
# BUILD A MOUSE SKIN scRNA-seq REFERENCE FOR RCTD
# OmicLearningHub · Module 3 case study · companion to xenium_tma_case_study.R
#
# Two sources (choose with SOURCE):
#   "GSE129218" – Joost et al. 2020, Cell Stem Cell: full-thickness dorsal
#                 skin (epidermis, hair follicle, dermis, immune, vascular),
#                 10x UMI counts, ~5,800 cells, 5- and 9-week-old mice.
#                 RECOMMENDED for full-thickness skin. You name clusters.
#   "TMS"       – Tabula Muris Senis, Skin (droplet / 10x). Ready-made
#                 cell_ontology_class labels, but check coverage first:
#                 if fibroblast / endothelial classes are missing, RCTD
#                 will force dermal cells into keratinocyte labels.
#
# Status: not yet test-run end to end. Each stage prints what it found and
# stops with a clear message if the data look different than expected.
# Re-runnable: downloads, loading and clustering are cached.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2)
})
options(timeout = 3600)          # large downloads
set.seed(1)

# ---------------- SETTINGS --------------------------------------------
SOURCE   <- "GSE129218"           # "GSE129218" (recommended) or "TMS"
REF_ROOT <- path.expand("~/refs")
OUT_RDS  <- file.path(REF_ROOT, paste0("mouse_skin_", SOURCE, "_rctd_ref.rds"))

# GSE129218 only – A) author annotation file, if you have one (else NULL)
ANN_FILE        <- NULL           # e.g. file.path(REF_ROOT, "GSE129218", "annotation.txt.gz")
ANN_BARCODE_COL <- NULL           # column with cell barcodes
ANN_CT_COL      <- NULL           # column with cell type names

# GSE129218 only – B) otherwise name clusters after the first run, using
#    ref_marker_dotplot.pdf and ref_cluster_markers.csv.
#    * clusters you leave out are dropped; the same name merges clusters
#    * keep labels COARSE (e.g. one "Immune", one "Vascular"): RCTD drops
#      every type with < 25 cells, and this dataset has only ~150-400
#      fibroblast / immune / vascular / neural-crest cells in total
#    * no comma after the last entry
CLUSTER_LABELS <- c(
  # "0" = "Basal_KC",
  # "1" = "Suprabasal_KC",
  # "2" = "Fibroblast"
)

# TMS only – keep these ages (NULL = all). Match your mice, e.g. c("1m", "3m")
TMS_AGES <- NULL

# Both – optional renaming / merging of final labels, old = new
# e.g. c("basal cell of epidermis" = "Basal_KC", "bulge keratinocyte" = "HairFollicle")
LABEL_MAP <- c()
# ----------------------------------------------------------------------

dir.create(REF_ROOT, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# SOURCE 1: GSE129218 (Joost et al. 2020)
# ======================================================================
if (SOURCE == "GSE129218") {
  suppressPackageStartupMessages(library(GEOquery))
  REF_DIR <- file.path(REF_ROOT, "GSE129218")
  dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)
  raw_rds <- file.path(REF_DIR, "ref_raw.rds")
  clu_rds <- file.path(REF_DIR, "ref_clustered.rds")

  # 1. DOWNLOAD (only if the folder is empty) ----------------------------
  if (length(list.files(REF_DIR, recursive = TRUE)) == 0) {
    message("Supplementary files on GEO:")
    print(getGEOSuppFiles("GSE129218", baseDir = REF_ROOT, fetch_files = FALSE))
    getGEOSuppFiles("GSE129218", baseDir = REF_ROOT)
  }
  # unpack every .tar into its own folder (all of them, not just the first)
  for (t in list.files(REF_DIR, pattern = "\\.tar$", full.names = TRUE)) {
    ex <- file.path(REF_DIR, "raw", sub("\\.tar$", "", basename(t)))
    if (!dir.exists(ex)) untar(t, exdir = ex)
  }
  message("Files in ", REF_DIR, ":"); print(list.files(REF_DIR, recursive = TRUE))

  # 2. LOAD 10x MATRICES -------------------------------------------------
  if (!file.exists(raw_rds)) {
    mtx <- list.files(REF_DIR, pattern = "matrix\\.mtx", recursive = TRUE, full.names = TRUE)
    if (length(mtx) == 0)
      stop("No matrix.mtx files found – the data are in another format (e.g. .loom / .h5 / .csv). ",
           "Copy the file list printed above and adapt this loading step.")
    objs <- lapply(mtx, function(m) {
      pre <- sub("matrix\\.mtx.*$", "", basename(m)); d <- dirname(m)
      esc <- gsub("([.])", "\\\\\\1", pre)
      bc  <- list.files(d, paste0("^", esc, ".*barcodes"), full.names = TRUE)
      ft  <- list.files(d, paste0("^", esc, ".*(genes|features)"), full.names = TRUE)
      if (length(bc) != 1 || length(ft) != 1) stop("Could not pair barcodes/features for ", m)
      proj <- if (nchar(pre) > 0) sub("_$", "", pre) else basename(d)
      CreateSeuratObject(ReadMtx(mtx = m, cells = bc, features = ft), project = proj)
    })
    ref <- if (length(objs) > 1) {
      JoinLayers(merge(objs[[1]], y = objs[-1],
                       add.cell.ids = sapply(objs, function(o) as.character(o$orig.ident[1]))))
    } else objs[[1]]
    rm(objs); gc()
    saveRDS(ref, raw_rds)
  } else ref <- readRDS(raw_rds)
  message("Loaded: ", ncol(ref), " cells x ", nrow(ref), " genes; samples: ",
          paste(unique(ref$orig.ident), collapse = ", "))

  # 3. QC ------------------------------------------------------------------
  ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^mt-")
  print(summary(ref@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")]))
  ref <- subset(ref, nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 15)
  message("After QC: ", ncol(ref), " cells")

  # 4. CELL TYPE LABELS --------------------------------------------------
  if (!is.null(ANN_FILE)) {
    # ---- A) labels from an annotation file
    ann <- fread(ANN_FILE, data.table = FALSE)
    if (is.null(ANN_BARCODE_COL) || is.null(ANN_CT_COL)) {
      print(head(ann)); print(head(colnames(ref)))
      stop("Set ANN_BARCODE_COL and ANN_CT_COL from the columns printed above, then re-run.")
    }
    k <- match(colnames(ref), ann[[ANN_BARCODE_COL]])
    if (mean(!is.na(k)) < 0.5) {                       # retry without sample prefix / -1 suffix
      strip <- function(x) sub("-1$", "", sub("^.*_", "", x))
      k <- match(strip(colnames(ref)), strip(ann[[ANN_BARCODE_COL]]))
    }
    message("Annotation matched for ", round(100 * mean(!is.na(k))), "% of cells")
    if (mean(!is.na(k)) < 0.5) stop("Barcodes don't match – compare head(ann) and head(colnames(ref)).")
    ref$celltype <- ann[[ANN_CT_COL]][k]
  } else {
    # ---- B) cluster the reference and name clusters by markers
    if (!file.exists(clu_rds)) {
      ref <- NormalizeData(ref) |> FindVariableFeatures(nfeatures = 3000) |> ScaleData() |> RunPCA(npcs = 30)
      ref <- FindNeighbors(ref, dims = 1:30) |> FindClusters(resolution = 0.5) |> RunUMAP(dims = 1:30)
      saveRDS(ref, clu_rds)
    } else ref <- readRDS(clu_rds)

    markers <- c("Krt14","Krt5","Krt10","Lor","Krt79","Krt17","Sox9","Lgr5","Krt25","Krt71",
                 "Scd1","Col1a1","Pdgfra","Dcn","Crabp1","Corin","Pecam1","Cdh5","Lyve1",
                 "Rgs5","Acta2","Ptprc","Cd3e","Trdc","Adgre1","Lyz2","Cd207","S100a8",
                 "Cpa3","Pmel","Dct","Mpz","Plp1","Adipoq","Acta1","Hba-a1")
    markers <- intersect(markers, rownames(ref))
    ggsave(file.path(REF_DIR, "ref_marker_dotplot.pdf"),
           DotPlot(ref, features = markers) + RotatedAxis(), width = 16, height = 8)
    ggsave(file.path(REF_DIR, "ref_umap.pdf"),
           DimPlot(ref, label = TRUE) | DimPlot(ref, group.by = "orig.ident"), width = 16, height = 7)
    if (!file.exists(file.path(REF_DIR, "ref_cluster_markers.csv"))) {
      top <- FindAllMarkers(ref, only.pos = TRUE, logfc.threshold = 0.5,
                            min.pct = 0.2, max.cells.per.ident = 300)
      write.csv(top, file.path(REF_DIR, "ref_cluster_markers.csv"), row.names = FALSE)
    }
    print(table(cluster = ref$seurat_clusters))

    if (length(CLUSTER_LABELS) == 0)
      stop("Clustering done. Open ", REF_DIR, "/ref_marker_dotplot.pdf and ref_cluster_markers.csv, ",
           "fill CLUSTER_LABELS at the top, and re-run.")
    ref$celltype <- unname(CLUSTER_LABELS[as.character(ref$seurat_clusters)])
  }
}

# ======================================================================
# SOURCE 2: Tabula Muris Senis, Skin (droplet / 10x)
# ======================================================================
if (SOURCE == "TMS") {
  suppressPackageStartupMessages({ library(TabulaMurisSenisData); library(SingleCellExperiment) })
  TabulaMurisSenisDroplet(tissues = "Skin", infoOnly = TRUE)        # prints download size
  sce <- TabulaMurisSenisDroplet(tissues = "Skin")$Skin             # raw counts (processedCounts = FALSE)
  print(colnames(colData(sce)))
  if (!"cell_ontology_class" %in% colnames(colData(sce)))
    stop("No 'cell_ontology_class' column – pick the label column from the list above.")
  print(table(celltype = sce$cell_ontology_class, age = sce$age))
  # FACS (Smart-seq2) data hold read counts, not UMIs; droplet suits RCTD better.

  if (!is.null(TMS_AGES)) sce <- sce[, sce$age %in% TMS_AGES]
  cnt <- as(counts(sce), "dgCMatrix")
  ref <- CreateSeuratObject(cnt, meta.data = as.data.frame(colData(sce)))
  ref$celltype <- as.character(ref$cell_ontology_class)
  rm(sce, cnt); gc()
}

if (!SOURCE %in% c("GSE129218", "TMS")) stop("SOURCE must be \"GSE129218\" or \"TMS\".")

# ======================================================================
# FINAL CHECKS & SAVE (both sources)
# ======================================================================
cnt_x <- LayerData(ref, "counts")@x
if (!all(cnt_x == round(cnt_x)))
  stop("Counts are not whole numbers – RCTD needs raw integer counts.")

if (length(LABEL_MAP) > 0) {
  hit <- ref$celltype %in% names(LABEL_MAP)
  ref$celltype[hit] <- unname(LABEL_MAP[ref$celltype[hit]])
}
ref <- subset(ref, cells = colnames(ref)[!is.na(ref$celltype) & ref$celltype != ""])
ref$celltype <- gsub("[/ ,]+", "_", ref$celltype)          # RCTD dislikes "/" and spaces

ct_n <- sort(table(ref$celltype), decreasing = TRUE); print(ct_n)
if (any(ct_n < 25))
  message("RCTD will drop these types (< 25 cells) – merge them into a coarser label: ",
          paste(names(ct_n)[ct_n < 25], collapse = ", "))
if (!any(grepl("fibro", names(ct_n), ignore.case = TRUE)) ||
    !any(grepl("endo|vasc", names(ct_n), ignore.case = TRUE)))
  warning("No fibroblast and/or endothelial label – this reference may not cover the dermis. ",
          "RCTD will force those cells into other types.")

ref <- DietSeurat(ref, layers = "counts")                 # RCTD needs counts only
saveRDS(ref, OUT_RDS)
message("Reference saved: ", OUT_RDS,
        "\nIn xenium_tma_case_study.R set:\n  RUN_RCTD <- TRUE\n  REF_RDS <- \"", OUT_RDS,
        "\"\n  REF_CT_COL <- \"celltype\"")
