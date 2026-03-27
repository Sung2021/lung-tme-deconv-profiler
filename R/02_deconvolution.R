# =============================================================================
# 02_deconvolution.R  [WRAPPER — backwards-compatibility shim]
#
# This file is kept so that direct calls to  Rscript R/02_deconvolution.R
# continue to work. The actual logic lives in two sub-scripts:
#   02a_signature_matrix.R  — load, validate, and prepare the signature
#   02b_deconvolution.R     — run QP + NNLS, compare, save results
#
# The Python orchestrator (python/run_pipeline.py) runs 02a and 02b directly.
# =============================================================================

source("R/02a_signature_matrix.R")
source("R/02b_deconvolution.R")


# --- Step 1: Load configuration ---
config <- yaml::read_yaml("config.yaml")

# --- Step 2: Load normalized expression and signature matrix ---
# !! SIGNATURE MATRIX NOTE !!
# The signature matrix (sig_matrix.csv) defines which genes distinguish each
# cell type (Cancer, Immune, Stroma). Its quality directly determines the
# accuracy of all downstream proportions and scores.
#
# Recommended sources:
#   - LM22    : 22-immune-cell signature, best validated for CIBERSORT-style QP
#               (Newman et al. 2015, Nature Methods)
#   - EPIC     : Cancer/Immune/Stroma 3-class, designed for solid tumours
#               (Racle et al. 2017, eLife) — closest match to this project
#   - MCP-counter : marker-gene based, no sum-to-one constraint
#
# Format requirements for sig_matrix.csv:
#   - Rows    : gene symbols (must match rownames of the expression matrix)
#   - Columns : cell type names (e.g., Cancer, Immune, Stroma)
#   - Values  : log2 TPM (must match the normalization used in 01_preprocessing.R)
#               If signature is in a different unit, re-scale before saving to
#               data/reference/sig_matrix.csv.
#
# CRITICAL: expression matrix and signature matrix MUST use the same:
#   (a) gene ID type  (HGNC symbol vs Ensembl ID)
#   (b) normalization (both log2 TPM here)
expr <- read.csv(config$paths$normalized_expr, row.names = 1, check.names = FALSE)
sig  <- read.csv(config$paths$signature_matrix, row.names = 1, check.names = FALSE)
cat("Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("Signature matrix :", nrow(sig),  "genes x", ncol(sig),  "cell types:", paste(colnames(sig), collapse = ", "), "\n")

# --- Step 3: Align genes and validate signature quality ---
common_genes <- intersect(rownames(expr), rownames(sig))
coverage_pct <- round(length(common_genes) / nrow(sig) * 100, 1)
cat(sprintf("Common genes: %d / %d signature genes  (%.1f%% coverage)\n",
            length(common_genes), nrow(sig), coverage_pct))

# Warn if coverage is too low to trust deconvolution results
if (coverage_pct < 50) {
  warning(sprintf(
    "Only %.1f%% of signature genes found in expression matrix. ",
    "Check that gene ID types and normalization methods match between ",
    "the expression matrix (01_preprocessing.R) and sig_matrix.csv.",
    coverage_pct
  ))
}

expr_aligned <- as.matrix(expr[common_genes, ])
sig_aligned  <- as.matrix(sig[common_genes, ])

# Check condition number of the aligned signature matrix
# A high condition number (> 1000) indicates collinearity between cell types,
# which makes the QP solution numerically unstable.
kappa_val <- kappa(sig_aligned, exact = TRUE)
cat(sprintf("Signature matrix condition number: %.1f %s\n",
            kappa_val,
            ifelse(kappa_val > 1000, "(WARNING: high collinearity — results may be unstable)",
                   "(OK)")))

# --- Step 4: QP deconvolution per sample ---
n_samples   <- ncol(expr_aligned)
n_celltypes <- ncol(sig_aligned)

# QP setup: minimize ||Sig * w - y||^2 subject to w >= 0, sum(w) = 1
Dmat <- t(sig_aligned) %*% sig_aligned
Amat <- cbind(rep(1, n_celltypes), diag(n_celltypes))    # equality + non-negativity
bvec <- c(1, rep(0, n_celltypes))                          # sum=1, each >= 0
meq  <- 1                                                  # first constraint is equality

proportions <- matrix(NA, nrow = n_samples, ncol = n_celltypes)
rownames(proportions) <- colnames(expr_aligned)
colnames(proportions) <- colnames(sig_aligned)

for (i in seq_len(n_samples)) {
  dvec <- t(sig_aligned) %*% expr_aligned[, i]
  sol  <- solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
  proportions[i, ] <- sol$solution
}

proportions_df <- as.data.frame(proportions)
proportions_df$SampleID <- rownames(proportions_df)
cat("Deconvolution complete for", n_samples, "samples.\n")

# --- Step 5: Save cell proportions ---
dir.create(dirname(config$paths$cell_proportions), recursive = TRUE, showWarnings = FALSE)
write.csv(proportions_df, config$paths$cell_proportions, row.names = FALSE)
cat("Cell proportions saved to:", config$paths$cell_proportions, "\n")

# --- Step 6: Plot stacked bar chart ---
plot_df <- tidyr::pivot_longer(
  proportions_df,
  cols      = -SampleID,
  names_to  = "CellType",
  values_to = "Proportion"
)

p <- ggplot(plot_df, aes(x = SampleID, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(title = "Cell Type Composition per Sample", x = "Sample", y = "Proportion")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/cell_composition_barplot.png", p, width = 12, height = 6, dpi = 300)
cat("Stacked bar chart saved to: results/figures/cell_composition_barplot.png\n")
