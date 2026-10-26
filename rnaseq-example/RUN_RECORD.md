# Run record · Module 1 (bulk RNA-seq)

- **Date:** 2026-10-26 22:36 UTC
- **Result:** 5 of 5 checks match the published workflow (Law et al. 2016)
- **Runtime:** 2:18.16 (m:ss) wall clock · **peak memory:** 1.76 GB
- **Machine:** Ubuntu 26.04.1 LTS, 4 CPUs, 15 GB RAM (GitHub Actions)
- **Commit:** `c19cc18` · **R:** 4.5.3 (2026-03-11)
- **Packages:** limma 3.66.0, edgeR 4.8.2, fgsea 1.36.2, clusterProfiler 4.18.4, msigdbr 26.1.1, org.Mm.eg.db 3.22.0, data.table 1.18.6.1
- **Environment:** [`environment_lock.yml`](expected/environment_lock.yml) (exact conda versions)

## Checks against the paper

| Check | Paper | This run | Status |
|---|---|---|---|
| genes in raw count table | 27179 | 27179 | match |
| genes kept by filterByExpr | 16624 | 16624 | match |
| TMM normalisation factors | 0.894 1.025 1.046 1.046 1.016 0.922 0.996 1.086 0.984 | 0.894 1.025 1.046 1.046 1.016 0.922 0.996 1.086 0.984 | match |
| eBayes DE counts (Down/NotSig/Up x 3 contrasts) | 4648 7113 4863 4927 7026 4671 3135 10972 2517 | 4646 7118 4860 4936 7008 4680 3141 10953 2530 | match |
| treat \|log2FC\|>1 DE counts | 1632 12976 2016 1777 12790 2057 224 16210 190 | 1633 12977 2014 1777 12793 2054 223 16211 190 | match |

DE counts differ by at most 19 genes (< 0.3 %) from the paper because limma/edgeR and the annotation have been updated since 2016; genes kept, TMM factors and the overall pattern are identical.

## Top Hallmark gene sets, Basal vs LP (GSEA, all genes ranked by t)

| Gene set | NES | adj. P |
|---|---|---|
| HALLMARK_E2F_TARGETS | 2.73 | 2.8e-22 |
| HALLMARK_G2M_CHECKPOINT | 2.67 | 1.3e-20 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 2.65 | 6.1e-20 |
| HALLMARK_ESTROGEN_RESPONSE_EARLY | -2.09 | 1.4e-08 |
| HALLMARK_MITOTIC_SPINDLE | 1.90 | 6.5e-06 |
| HALLMARK_MYC_TARGETS_V1 | 1.82 | 3.6e-05 |
