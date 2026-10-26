# OmicLearningHub

### 🌐 **[Open the website → githubrudramani.github.io/OmicLearningHub](https://githubrudramani.github.io/OmicLearningHub/)**

Free, hands-on bioinformatics tutorials. Each module explains how the technology developed, compares the main platforms, explains the key method with an interactive demo, and walks through a complete analysis with code you can copy on public data. Each module also says how far its code has been tested.

## Learning path

| # | Module | Covers | Status |
|---|---|---|---|
| 1 | [Bulk RNA-seq](https://githubrudramani.github.io/OmicLearningHub/rnaseq.html) | limma-voom + TMM on GSE63310, camera, MSigDB GSEA (all genes ranked by t) and ORA · [GSEA script](scripts/rnaseq_gsea_msigdb.R) · [full reproducible run](rnaseq-example/) | **Run on public data** (26 Oct 2026, 5 of 5 checks match the paper, [run record](rnaseq-example/RUN_RECORD.md)) |
| 2 | [Single-cell RNA-seq](https://githubrudramani.github.io/OmicLearningHub/scrnaseq.html) | LogNormalize, Harmony, ScType + SingleR, CellTypist alternative, paired pseudobulk DE · [CellTypist script](scripts/scrnaseq_celltypist.py) | Tested on simulated data (pseudobulk, CellTypist); full pipeline not yet run |
| 3 | [Spatial (Xenium)](https://githubrudramani.github.io/OmicLearningHub/xenium.html) | Platforms, Seurat Xenium pipeline, Banksy, niches; [multi-slide TMA case study](https://githubrudramani.github.io/OmicLearningHub/xenium.html#case-study) · [script](scripts/xenium_tma_case_study.R) · [RCTD reference builder](scripts/build_skin_rctd_reference.R) | Tested on simulated data (sketch projection, case-study checks); full pipeline not yet run |

Coming next: Visium HD · DNA-seq / WGS · ATAC-seq

**Validation labels:** *Run on public data*: the whole pipeline was executed on the named public dataset, with outputs and versions published. *Tested on simulated data*: code run on simulated data with known answers. *Not yet run*: written from official documentation; adapt and check it.

## Report a problem · cite

Open an [issue](https://github.com/githubrudramani/OmicLearningHub/issues/new) naming the module and step.
Cite as: githubrudramani (2026). *OmicLearningHub: hands-on bioinformatics tutorials.* https://githubrudramani.github.io/OmicLearningHub/, and also cite the methods you use.

## How the site is built

Plain HTML, CSS and JavaScript with no build step. Every push to `main` is published by `.github/workflows/pages.yml`. The bulk RNA-seq module is re-run end to end by `.github/workflows/rnaseq-example.yml` whenever its code changes; results are pushed to the [`rnaseq-run`](https://github.com/githubrudramani/OmicLearningHub/tree/rnaseq-run) branch.
