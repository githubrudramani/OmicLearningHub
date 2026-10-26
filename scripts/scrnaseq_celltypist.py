#!/usr/bin/env python3
# ======================================================================
# CELLTYPIST ANNOTATION FOR A SEURAT OBJECT
# OmicLearningHub · Module 3 · https://githubrudramani.github.io/OmicLearningHub/scrnaseq.html
#
# Input: the folder written by the R export snippet (step 13 on the page)
#   counts.mtx   genes x cells raw counts (Matrix::writeMM)
#   genes.txt    gene symbols, one per line
#   cells.txt    cell barcodes, one per line
#   clusters.csv cell,cluster  (Harmony clusters, used for majority voting)
# Output (same folder): celltypist_<model>.csv with, per cell,
#   predicted_labels, majority_voting, conf_score
#
# Usage:
#   pip install celltypist scanpy
#   (if celltypist.train() fails with "unexpected keyword argument
#    'multi_class'", install "scikit-learn<1.7": CellTypist 1.7.1 predates it)
#   python scrnaseq_celltypist.py --indir celltypist \
#          --models Immune_All_Low.pkl Immune_All_High.pkl
#   (a path to your own .pkl model also works, e.g. one made with celltypist.train)
# ======================================================================
import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
import scipy.io
import celltypist
from celltypist import models

ap = argparse.ArgumentParser(description="Annotate exported Seurat counts with CellTypist")
ap.add_argument("--indir", default="celltypist", help="folder written by the R export step")
ap.add_argument("--models", nargs="+", default=["Immune_All_Low.pkl", "Immune_All_High.pkl"],
                help="built-in model names or paths to .pkl files")
ap.add_argument("--no-majority-voting", action="store_true",
                help="skip majority voting (labels per cell only)")
args = ap.parse_args()
d = Path(args.indir)

# 1. Load counts as cells x genes ------------------------------------------
X = scipy.io.mmread(d / "counts.mtx").T.tocsr()
genes = (d / "genes.txt").read_text().split("\n")
cells = (d / "cells.txt").read_text().split("\n")
genes = [g for g in genes if g]; cells = [c for c in cells if c]
adata = sc.AnnData(X)
adata.var_names = genes
adata.obs_names = cells
adata.var_names_make_unique()
clusters = pd.read_csv(d / "clusters.csv", index_col=0)
adata.obs["cluster"] = clusters.loc[adata.obs_names, "cluster"].astype(str).values
print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes, {adata.obs['cluster'].nunique()} clusters")

# 2. CellTypist needs log1p(counts per 10,000) – the same as Seurat LogNormalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 3. Get built-in models (skips anything given as a file path) ------------
builtin = [m for m in args.models if not Path(m).exists()]
if builtin:
    models.download_models(model=builtin)    # downloads once to ~/.celltypist

# 4. Annotate --------------------------------------------------------------
for m in args.models:
    # local files must be passed with a path: CellTypist treats a bare name
    # without "/" as a built-in model and tries to download all models
    model = str(Path(m).resolve()) if Path(m).exists() else m
    pred = celltypist.annotate(adata, model=model,
                               majority_voting=not args.no_majority_voting,
                               over_clustering="cluster")
    res = pred.to_adata().obs
    cols = [c for c in ["predicted_labels", "majority_voting", "conf_score"] if c in res.columns]
    out = d / f"celltypist_{Path(m).stem}.csv"
    res[cols].to_csv(out, index_label="cell")
    label = "majority_voting" if "majority_voting" in cols else "predicted_labels"
    print(f"\n{Path(m).stem}: wrote {out}")
    print(res[label].value_counts().head(15).to_string())
    print(f"median confidence: {res['conf_score'].median():.2f}; "
          f"cells with conf_score < 0.5: {(res['conf_score'] < 0.5).mean():.1%}")
