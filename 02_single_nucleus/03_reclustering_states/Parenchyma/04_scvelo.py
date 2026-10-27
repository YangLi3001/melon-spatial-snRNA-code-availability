from pathlib import Path

import anndata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvelo as scv

# Inputs: 03_prepare_velocity_h5ad.R output and the study's merged loom file.
sce = anndata.read_h5ad("parenchyma.h5ad")
ldat = sc.read("velocity.loom", cache=True, sparse=True)

# Remove the library prefixes used in the original loom barcodes.
prefixes = ["sorted_velocity_", "YMCUT:", "GZLTP:", "KLAHY:", "WJXCF:"]
barcodes = ldat.obs_names.astype(str).tolist()
for prefix in prefixes:
    barcodes = [barcode.replace(prefix, "") for barcode in barcodes]
ldat.obs_names = barcodes
ldat.var_names_make_unique()
if not ldat.obs_names.is_unique:
    raise ValueError("Duplicate loom barcodes after prefix removal")
if not set(sce.obs_names).intersection(ldat.obs_names):
    raise ValueError("No matching barcodes between the h5ad and loom files")
adata = scv.utils.merge(sce, ldat)

scv.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
scv.pp.moments(adata, n_neighbors=30, n_pcs=30)
scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata, n_jobs=40)
scv.tl.velocity_confidence(adata)
scv.tl.velocity_embedding(adata, basis="umap")

out_dir = Path("results/velocity")
out_dir.mkdir(parents=True, exist_ok=True)
adata.write(out_dir / "velocity.h5ad", compression="gzip")

# Retain finite coordinates for the velocity projection.
finite = np.isfinite(adata.obsm["X_umap"]).all(axis=1)
finite &= np.isfinite(adata.obsm["velocity_umap"]).all(axis=1)
adata_plot = adata[finite].copy()
scv.pl.velocity_embedding_stream(adata_plot, basis="umap", color="celltype", density=2, linewidth=1.2, show=False)
plt.savefig(out_dir / "velocity_stream_umap.png", dpi=300, bbox_inches="tight")
plt.close("all")
