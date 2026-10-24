# OmicsRun Learn

Static learning hub for bioinformatics pipelines, served with GitHub Pages.

| Page | Module |
|---|---|
| `index.html` | Home |
| `xenium.html` | Module 1 – Xenium spatial transcriptomics (Seurat, Harmony, Banksy) |
| `rnaseq.html` | Module 2 – Bulk RNA-seq with limma-voom and TMM |
| `scrnaseq.html` | Module 3 – scRNA-seq with LogNormalize, Harmony, ScType/SingleR |

## Deploy
1. Push this folder to a GitHub repository (branch `main`).
2. In the repo: **Settings → Pages → Build and deployment → Source: GitHub Actions**.
3. The workflow in `.github/workflows/pages.yml` publishes the site to
   `https://<your-username>.github.io/<repo-name>/`.

No build step is needed: plain HTML, CSS and JavaScript.
