# OmicLearningHub

### 🌐 **[Open the website → githubrudramani.github.io/OmicLearningHub](https://githubrudramani.github.io/OmicLearningHub/)**

Hands-on bioinformatics learning hub. Each module explains how the technology developed, compares the main platforms, explains the key method with an interactive demo, and walks through a complete analysis with code you can copy on public data.

| Module | Topic | Open |
|---|---|---|
| 1 · Spatial | Xenium spatial transcriptomics: platform comparison, Seurat, Harmony, Banksy, niches | [xenium.html](https://githubrudramani.github.io/OmicLearningHub/xenium.html) |
| 2 · Bulk RNA-seq | limma-voom with TMM normalisation on GSE63310 | [rnaseq.html](https://githubrudramani.github.io/OmicLearningHub/rnaseq.html) |
| 3 · Single-cell | scRNA-seq with LogNormalize, Harmony integration, ScType + SingleR annotation | [scrnaseq.html](https://githubrudramani.github.io/OmicLearningHub/scrnaseq.html) |

Coming next: Visium HD · DNA-seq / WGS · ATAC-seq

## How the site is built
Plain HTML, CSS and JavaScript with no build step. Every push to `main` is published by the GitHub Actions workflow in `.github/workflows/pages.yml`.
