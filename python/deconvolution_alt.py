#!/usr/bin/env python3
"""
deconvolution_alt.py
--------------------
Alternative TME deconvolution methods using Python (scipy / sklearn).

Two methods are implemented and compared against each other:

  Method 1 — NNLS (Non-Negative Least Squares, supervised, reference-based)
      Solves per sample:  min ||S·w - y||²  s.t.  w >= 0
      then normalizes each solution to sum to 1.
      Functionally equivalent to the QP approach in 02_deconvolution.R
      but uses scipy's faster NNLS solver. No sum-to-one constraint is
      enforced during optimization (post-hoc normalization instead).

  Method 2 — NMF (Non-negative Matrix Factorization, unsupervised, reference-free)
      Factorizes the full expression matrix X ≈ W × H without a reference
      signature. Components are labeled by matching them to reference cell
      types via Pearson correlation (requires signature matrix for labeling
      only; not for the factorization itself).

Inputs (paths from config.yaml):
    - normalized expression matrix (genes × samples)
    - reference signature matrix   (genes × cell_types)

Outputs:
    - data/processed/cell_proportions_nnls.csv
    - data/processed/cell_proportions_nmf.csv
    - results/figures/deconvolution_comparison.png

Requirements:
    pip install pyyaml numpy pandas scipy scikit-learn matplotlib
"""

import warnings
import yaml
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend — safe for headless execution
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import nnls
from sklearn.decomposition import NMF
from sklearn.preprocessing import normalize

warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Step 1: Load configuration and resolve paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent

with open(ROOT / "config.yaml") as f:
    config = yaml.safe_load(f)

expr_path = ROOT / config["paths"]["normalized_expr"]
sig_path  = ROOT / config["paths"]["signature_matrix"]

print("TME Deconvolution — Alternative Methods (NNLS + NMF)")
print(f"Expression : {expr_path}")
print(f"Signature  : {sig_path}")

# ---------------------------------------------------------------------------
# Step 2: Load and align matrices
# ---------------------------------------------------------------------------
expr = pd.read_csv(expr_path, index_col=0)   # genes × samples
sig  = pd.read_csv(sig_path,  index_col=0)   # genes × cell_types

common_genes = expr.index.intersection(sig.index)
print(f"Common genes used : {len(common_genes)}")

X = expr.loc[common_genes].values     # (G × N)  genes × samples
S = sig.loc[common_genes].values      # (G × K)  genes × cell_types

cell_types = sig.columns.tolist()
samples    = expr.columns.tolist()

N = len(samples)
K = len(cell_types)

# ---------------------------------------------------------------------------
# Step 3: Method 1 — NNLS (supervised, per-sample)
# ---------------------------------------------------------------------------
print(f"\n[Method 1] NNLS deconvolution ({N} samples × {K} cell types)...")

nnls_matrix = np.zeros((N, K))

for i in range(N):
    coef, _ = nnls(S, X[:, i])
    total = coef.sum()
    nnls_matrix[i] = coef / total if total > 0 else coef  # normalize to sum=1

nnls_df = pd.DataFrame(nnls_matrix, index=samples, columns=cell_types)
nnls_df.index.name = "SampleID"
nnls_df = nnls_df.reset_index()

out_nnls = ROOT / "data" / "processed" / "cell_proportions_nnls.csv"
out_nnls.parent.mkdir(parents=True, exist_ok=True)
nnls_df.to_csv(out_nnls, index=False)
print(f"  Saved: {out_nnls}")

# ---------------------------------------------------------------------------
# Step 4: Method 2 — NMF (unsupervised)
# ---------------------------------------------------------------------------
print(f"\n[Method 2] NMF deconvolution (n_components={K}, random_state=42)...")

nmf_model = NMF(
    n_components=K,
    init="nndsvda",   # deterministic init — better for sparse data
    random_state=42,
    max_iter=1000,
    tol=1e-4,
)

# X.T shape: samples × genes  →  W: samples × components
W = nmf_model.fit_transform(X.T)          # (N × K)
H = nmf_model.components_                 # (K × G)
W_norm = normalize(W, norm="l1")          # row-normalize so each sample sums to 1

# Label NMF components by matching to reference cell types (highest Pearson r)
component_to_cell = {}
used_cell_types = set()

for comp_idx in range(K):
    comp_vec = H[comp_idx]
    correlations = {
        ct: np.corrcoef(comp_vec, S[:, j])[0, 1]
        for j, ct in enumerate(cell_types)
        if ct not in used_cell_types
    }
    if correlations:
        best_ct = max(correlations, key=correlations.get)
        used_cell_types.add(best_ct)
    else:
        best_ct = f"Component_{comp_idx + 1}"
    component_to_cell[comp_idx] = best_ct

component_labels = [component_to_cell[i] for i in range(K)]
print(f"  NMF component label mapping: {dict(enumerate(component_labels))}")

nmf_df = pd.DataFrame(W_norm, index=samples, columns=component_labels)
nmf_df.index.name = "SampleID"
nmf_df = nmf_df.reset_index()

out_nmf = ROOT / "data" / "processed" / "cell_proportions_nmf.csv"
nmf_df.to_csv(out_nmf, index=False)
print(f"  Saved: {out_nmf}")

# ---------------------------------------------------------------------------
# Step 5: Comparison — print mean proportions per method
# ---------------------------------------------------------------------------
print("\n--- NNLS  mean proportions ---")
print(nnls_df.drop(columns="SampleID").mean().round(4).to_string())

print("\n--- NMF   mean proportions ---")
print(nmf_df.drop(columns="SampleID").mean().round(4).to_string())

# ---------------------------------------------------------------------------
# Step 6: Stacked bar comparison (first 20 samples)
# ---------------------------------------------------------------------------
n_plot = min(20, N)
plot_samples = samples[:n_plot]

nnls_plot = nnls_df.set_index("SampleID").loc[plot_samples]
nmf_plot  = nmf_df.set_index("SampleID").loc[plot_samples]

fig, axes = plt.subplots(1, 2, figsize=(18, 6))
colors = plt.cm.Set2.colors  # type: ignore[attr-defined]

for ax, df, title in zip(axes, [nnls_plot, nmf_plot], ["NNLS (Supervised)", "NMF (Unsupervised)"]):
    df.plot(kind="bar", stacked=True, ax=ax, color=list(colors[:K]))
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Proportion")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=7)
    ax.legend(title="Cell Type", bbox_to_anchor=(1.01, 1.0), loc="upper left")
    ax.set_ylim(0, 1.05)

plt.suptitle(
    f"TME Deconvolution Comparison: NNLS vs NMF (first {n_plot} samples)",
    fontsize=14, y=1.01
)
plt.tight_layout()

out_fig = ROOT / "results" / "figures" / "deconvolution_comparison.png"
out_fig.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out_fig, dpi=300, bbox_inches="tight")
print(f"\nComparison plot saved: {out_fig}")
print("\nDone.")
