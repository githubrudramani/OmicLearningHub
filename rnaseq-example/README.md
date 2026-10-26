# Module 1 · reproducible bulk RNA-seq run

One command reproduces every result in [Module 1](https://githubrudramani.github.io/OmicLearningHub/rnaseq.html)
on the public GSE63310 dataset (Law et al. 2016): download → filtering → TMM → limma-voom → DE tables and plots →
camera → MSigDB GSEA / ORA. **Status: run on public data, 5 of 5 checks match the published workflow** (see `RUN_RECORD.md`).

## Run it

```bash
git clone https://github.com/githubrudramani/OmicLearningHub.git
cd OmicLearningHub
micromamba create -f rnaseq-example/environment.yml      # or: conda env create -f ...
micromamba activate omiclearninghub-rnaseq
Rscript rnaseq-example/run_rnaseq_limma.R output         # ~2–5 min, < 2 GB RAM
```

Without conda: install R ≥ 4.5, run `Rscript rnaseq-example/install.R`, then the same `Rscript` line.
The first run needs internet (GEO data ≈ 4 MB, MSigDB gene sets, org.Mm.eg.db).

## What you get (in `output/`)

| File | Content |
|---|---|
| `checks_vs_paper.csv` | Key numbers compared with the published workflow; every row should say `match` |
| `06_samples_libsize_normfactors.csv` | Library sizes and TMM factors |
| `09_decideTests_eBayes.csv`, `10_decideTests_treat_lfc1.csv` | DE counts per contrast |
| `10_DE_<contrast>_treat.csv` | Full DE tables (all genes) for each contrast |
| `11_camera_Hallmark_BasalvsLP.csv` | camera gene-set test |
| `gsea_results/` | GSEA (all genes ranked by t), ORA and restricted GSEA for Hallmark, Reactome, GO:BP |
| `plots/` | Density, MDS, voom trend, SA plot, MD plot, volcano, heatmap, barcode plot |
| `sessionInfo.txt`, `run_summary.txt` | Versions, runtime, check summary |

Reference copies of the small outputs and all plots from the last verified run are in [`expected/`](expected/).
Full outputs: branch [`rnaseq-run`](https://github.com/githubrudramani/OmicLearningHub/tree/rnaseq-run).

## Troubleshooting

| Problem | Fix |
|---|---|
| `Could not solve for environment specs` | Use the file unchanged; it pins only R and lets conda choose one consistent Bioconductor release. Mixing versions by hand (e.g. limma from one release, Glimma from another) fails. |
| `download.file` timeout / GEO error | GEO is sometimes slow. Re-run: the script skips files it already downloaded (`output/data/`). |
| msigdbr: `Downloading gene sets ... error` | Recent msigdbr downloads MSigDB from Zenodo on first use; you need internet once. Behind a proxy, set `https_proxy`. |
| Warning `cannot rename file ... Invalid cross-device link` (msigdbr) | Harmless: the gene sets still load; only the cache copy failed. |
| `calcNormFactors has been renamed to normLibSizes` | Informational message from edgeR 4; both names work. |
| A check says `DIFFERS` | Small shifts in DE counts are expected with newer package versions (tolerance ±100 genes). A large difference, or different genes kept / TMM factors, means the input or code differs: compare `sessionInfo.txt` with `expected/sessionInfo.txt`. |
| `there is no package called 'GO.db'` | clusterProfiler needs GO.db; with conda it comes with the environment, otherwise `BiocManager::install("GO.db")`. |

Report other problems as a [GitHub issue](https://github.com/githubrudramani/OmicLearningHub/issues/new?title=%5BModule%201%20reproducible%20run%5D%20).
