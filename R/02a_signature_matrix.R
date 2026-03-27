# =============================================================================
# 02a_signature_matrix.R
# Load, validate, and prepare the signature matrix for QP deconvolution
#
# WHY a separate script?
#   The signature matrix is the most critical input to deconvolution.
#   A poorly calibrated signature propagates directly into all cell proportion
#   estimates and downstream scores. Isolating this logic here allows the
#   matrix to be reviewed, swapped, or extended without touching the solver.
#
# DEFAULT SOURCE: EPIC R package (Racle et al. 2017, eLife)
#   EPIC::TRef contains 8 cell types calibrated to log2 TPM values:
#     Bcells, CAFs, CD4_Tcells, CD8_Tcells, Endothelial,
#     Macrophages, NKcells, otherCells
#   We consolidate to 3 coarse classes matching this project's scoring:
#     Cancer  <- otherCells
#     Immune  <- Bcells + CD4_Tcells + CD8_Tcells + Macrophages + NKcells  (mean)
#     Stroma  <- CAFs + Endothelial  (mean)
#
# ALTERNATIVE SOURCES (modify Step 2b):
#   - LM22 (Newman et al. 2015, Nat Methods): 22 immune subtypes, requires
#     registration at https://cibersortx.stanford.edu/
#   - Custom scRNA-seq: extract top marker genes per cell cluster from a
#     lung cancer atlas (e.g., HLCA — GSE136831) using Seurat::FindMarkers()
#     then construct pseudo-bulk mean expression per cluster as signature.
#
# FALLBACK: if EPIC is not installed, reads config$paths$signature_matrix
#
# Output: data/reference/sig_matrix_validated.csv  (genes x 3 cell types)
# =============================================================================

set.seed(42)

# --- Step 1: Load configuration ---
config <- yaml::read_yaml("config.yaml")

# --- Step 2: Load signature matrix ---
use_epic <- requireNamespace("EPIC", quietly = TRUE)

if (use_epic) {
  cat("[Signature] Source: EPIC package (TRef — Racle et al. 2017)\n")
  sig_raw <- as.data.frame(EPIC::TRef$refProfiles)
  cat("  Original EPIC cell types:", paste(colnames(sig_raw), collapse = ", "), "\n\n")

  # --- Step 3: Coarse labeling — 8 EPIC types → 3 classes ---
  # Averaging within each class produces a representative centroid profile.
  # This is appropriate because within-class profiles share strong correlation
  # (all immune cells are co-regulated) and averaging reduces noise.
  immune_cols <- intersect(colnames(sig_raw),
                           c("Bcells", "CD4_Tcells", "CD8_Tcells", "Macrophages", "NKcells"))
  stroma_cols <- intersect(colnames(sig_raw), c("CAFs", "Endothelial"))
  cancer_cols <- intersect(colnames(sig_raw), c("otherCells"))

  sig_coarse <- data.frame(
    Cancer = if (length(cancer_cols) > 0) rowMeans(sig_raw[, cancer_cols, drop = FALSE]) else stop("Cancer column not found in EPIC TRef"),
    Immune = if (length(immune_cols) > 0) rowMeans(sig_raw[, immune_cols, drop = FALSE]) else stop("No immune columns found in EPIC TRef"),
    Stroma = if (length(stroma_cols) > 0) rowMeans(sig_raw[, stroma_cols, drop = FALSE]) else stop("No stroma columns found in EPIC TRef"),
    row.names = rownames(sig_raw)
  )
  cat(sprintf(
    "  Coarse mapping:\n    Cancer (%d col): %s\n    Immune (%d cols): %s\n    Stroma (%d cols): %s\n\n",
    length(cancer_cols), paste(cancer_cols, collapse = "+"),
    length(immune_cols), paste(immune_cols, collapse = "+"),
    length(stroma_cols), paste(stroma_cols, collapse = "+")
  ))

} else {
  # --- Step 2b: Fallback — load CSV from config path ---
  cat("[Signature] EPIC not installed. Loading from:", config$paths$signature_matrix, "\n")
  cat("  For best results: install.packages('EPIC') and re-run.\n\n")
  sig_coarse <- read.csv(config$paths$signature_matrix, row.names = 1, check.names = FALSE)
}

cat("Signature dimensions:", nrow(sig_coarse), "genes x", ncol(sig_coarse), "cell types\n\n")

# --- Step 4: Quality validation ---
cat("=== Signature Matrix Quality Report ===\n\n")

# 4a: Per-column variance
# Each cell type column must have sufficient variance to serve as a discriminative
# feature. Zero variance means a column is constant — useless for regression.
col_vars <- apply(sig_coarse, 2, var, na.rm = TRUE)
cat("Per-cell-type variance (should all be > 0):\n")
print(round(col_vars, 5))
cat("\n")

zero_var_cols <- names(col_vars[col_vars < 1e-8])
if (length(zero_var_cols) > 0) {
  stop(paste("Zero-variance columns detected:", paste(zero_var_cols, collapse = ", "),
             "\nThese carry no discriminative signal. Check coarse labeling in Step 3."))
}

# 4b: Condition number
# Measures the numerical sensitivity of the signature matrix S.
# High kappa → small perturbations in Y produce large swings in the QP solution.
# Rule of thumb: < 100 excellent | 100–1000 acceptable | > 1000 unstable
kappa_val   <- kappa(as.matrix(sig_coarse), exact = TRUE)
kappa_label <- if      (kappa_val < 100)  "excellent"     \
               else if (kappa_val < 1000) "acceptable"    \
               else                       "WARNING: high — QP solutions may be unstable"
cat(sprintf("Condition number: %.1f  [%s]\n\n", kappa_val, kappa_label))

# 4c: Pairwise Pearson correlation between cell type profiles
# r > 0.85 between two types indicates they share so many marker genes that
# the solver cannot reliably distinguish them.
sig_cor <- cor(sig_coarse, use = "complete.obs")
cat("Pairwise Pearson correlation between cell type profiles:\n")
print(round(sig_cor, 3))
cat("\n")

high_cor_idx <- which(abs(sig_cor) > 0.85 & lower.tri(sig_cor), arr.ind = TRUE)
if (nrow(high_cor_idx) > 0) {
  pairs <- apply(high_cor_idx, 1, function(idx)
    paste(rownames(sig_cor)[idx[1]], "&", colnames(sig_cor)[idx[2]]))
  warning(paste("High correlation (r > 0.85) detected:", paste(pairs, collapse = "; "),
                "— deconvolution may not reliably separate these cell types."))
} else {
  cat("No high-correlation pairs detected (r <= 0.85 for all pairs). OK.\n\n")
}

# --- Step 5: Save validated signature matrix ---
out_path <- config$paths$sig_matrix_validated
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(sig_coarse, out_path, row.names = TRUE)
cat("Validated signature matrix saved to:", out_path, "\n")
