# =============================================================================
# 02b_deconvolution.R
# Run QP and NNLS deconvolution; compare results; save canonical output
#
# WHY two methods?
#   Both solve: min ||S*w - y||²  (least-squares fit of signature to expression)
#   but with different constraints on w (cell type proportions):
#
#   QP   [PRIMARY]   w_i >= 0  AND  Σw_i = 1  — enforced simultaneously
#   NNLS [REFERENCE] w_i >= 0  only            — sum-to-1 applied post-hoc
#
#   QP is the primary method because:
#   (1) sum-to-one is a biological hard constraint, not an approximation
#   (2) enforcing it during optimization avoids post-hoc normalization errors
#   (3) it is the mathematical basis of CIBERSORT (which uses SVR, a related
#       constrained solver but requires a commercial licence)
#
#   NNLS serves as a sanity check: if QP and NNLS agree (r > 0.95 per cell
#   type), results are numerically robust. Large discrepancies signal
#   collinearity or instability in the signature — revisit 02a.
#
# Outputs:
#   data/processed/proportions_qp.csv           (QP results)
#   data/processed/proportions_nnls.csv         (NNLS results)
#   data/processed/cell_proportions.csv         (QP — canonical, used downstream)
#   results/figures/cell_composition_barplot.png
#   results/figures/method_comparison.png
# =============================================================================

set.seed(42)
library(quadprog)
library(nnls)
library(ggplot2)

# --- Step 1: Load configuration ---
config <- yaml::read_yaml("config.yaml")

# --- Step 2: Load normalized expression and validated signature ---
expr <- read.csv(config$paths$normalized_expr,     row.names = 1, check.names = FALSE)
sig  <- read.csv(config$paths$sig_matrix_validated, row.names = 1, check.names = FALSE)
cat("Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("Signature matrix :", nrow(sig),  "genes x", ncol(sig),
    "cell types:", paste(colnames(sig), collapse = ", "), "\n")

# --- Step 3: Align genes and check coverage ---
common_genes <- intersect(rownames(expr), rownames(sig))
coverage_pct <- round(length(common_genes) / nrow(sig) * 100, 1)
cat(sprintf("Common genes: %d / %d signature genes (%.1f%% coverage)\n\n",
            length(common_genes), nrow(sig), coverage_pct))
if (coverage_pct < 50) {
  warning(sprintf(
    "Only %.1f%% coverage — check gene ID type and normalization units match between scripts.",
    coverage_pct
  ))
}

X          <- as.matrix(expr[common_genes, ])   # G x N: genes x samples
S          <- as.matrix(sig[common_genes, ])    # G x K: genes x cell types
N          <- ncol(X)
K          <- ncol(S)
cell_types <- colnames(S)

# --- Step 4: Method A — QP (Primary) ---
# Solve per sample: argmin_w ||S*w - y||²  s.t.  w_i >= 0, Σw_i = 1
# quadprog::solve.QP notation:
#   Dmat = S'S,  dvec = S'y,  Amat = [1|I],  bvec = [1|0...0],  meq = 1
cat("[Method A] QP deconvolution...\n")

Dmat <- t(S) %*% S
Amat <- cbind(rep(1, K), diag(K))   # first column: equality (sum=1); rest: non-negativity
bvec <- c(1, rep(0, K))
props_qp <- matrix(NA, nrow = N, ncol = K,
                   dimnames = list(colnames(X), cell_types))

for (i in seq_len(N)) {
  dvec         <- t(S) %*% X[, i]
  props_qp[i, ] <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1L)$solution
}
cat("  Done. Mean proportions (QP):\n")
print(round(colMeans(props_qp), 4))

# --- Step 5: Method B — NNLS (Reference) ---
# Solve per sample: argmin_w ||S*w - y||²  s.t.  w_i >= 0  (no equality constraint)
# Post-hoc row-normalize to sum = 1 for fair comparison.
cat("\n[Method B] NNLS deconvolution...\n")

props_nnls <- matrix(NA, nrow = N, ncol = K,
                     dimnames = list(colnames(X), cell_types))
for (i in seq_len(N)) {
  coef           <- nnls::nnls(S, X[, i])$x
  total          <- sum(coef)
  props_nnls[i, ] <- if (total > 0) coef / total else coef   # normalize to sum=1
}
cat("  Done. Mean proportions (NNLS):\n")
print(round(colMeans(props_nnls), 4))

# --- Step 6: Method comparison ---
cat("\n=== Method Comparison Summary ===\n")
cat("Pearson r (QP vs NNLS) per cell type:\n")
r_vals <- sapply(cell_types, function(ct)
  cor(props_qp[, ct], props_nnls[, ct], method = "pearson"))
print(round(r_vals, 4))

mean_abs_diff <- colMeans(abs(props_qp - props_nnls))
cat("\nMean absolute difference (QP - NNLS) per cell type:\n")
print(round(mean_abs_diff, 4))

if (all(r_vals > 0.95)) {
  cat("\n[PASS] Both methods agree closely (r > 0.95 for all cell types).\n")
  cat("       QP results are numerically robust and used as canonical output.\n\n")
} else {
  low_types <- names(r_vals[r_vals <= 0.95])
  warning(paste("Low agreement (r <= 0.95) for:", paste(low_types, collapse = ", "),
                "— re-check signature matrix (run 02a) or inspect gene coverage."))
}

# --- Step 7: Save all results ---
dir.create(dirname(config$paths$cell_proportions), recursive = TRUE, showWarnings = FALSE)

qp_df   <- as.data.frame(props_qp);   qp_df$SampleID   <- rownames(qp_df)
nnls_df <- as.data.frame(props_nnls); nnls_df$SampleID  <- rownames(nnls_df)

write.csv(qp_df,   config$paths$proportions_qp,   row.names = FALSE)
write.csv(nnls_df, config$paths$proportions_nnls,  row.names = FALSE)
write.csv(qp_df,   config$paths$cell_proportions,  row.names = FALSE)  # canonical
cat("Saved:\n",
    "  QP       :", config$paths$proportions_qp, "\n",
    "  NNLS     :", config$paths$proportions_nnls, "\n",
    "  Canonical:", config$paths$cell_proportions, "\n\n")

# --- Step 8: Plot A — stacked bar chart of QP proportions ---
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

bar_df <- tidyr::pivot_longer(qp_df, cols = -SampleID,
                              names_to = "CellType", values_to = "Proportion")
p_bar <- ggplot(bar_df, aes(x = SampleID, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(Cancer = "#E74C3C", Immune = "#2ECC71", Stroma = "#3498DB")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(title = "Cell Type Composition per Sample (QP — primary)",
       x = "Sample", y = "Proportion")
ggsave("results/figures/cell_composition_barplot.png", p_bar, width = 12, height = 6, dpi = 300)
cat("Stacked bar saved: results/figures/cell_composition_barplot.png\n")

# --- Step 9: Plot B — QP vs NNLS scatter per cell type ---
comp_df <- do.call(rbind, lapply(cell_types, function(ct) {
  data.frame(CellType = ct, QP = props_qp[, ct], NNLS = props_nnls[, ct])
}))
r_labels <- data.frame(CellType = cell_types,
                       label    = sprintf("r = %.3f", r_vals))

p_comp <- ggplot(comp_df, aes(x = QP, y = NNLS)) +
  geom_point(alpha = 0.5, size = 1.5, color = "#555555") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_text(data = r_labels, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.15, vjust = 1.8,
            size = 3.5, color = "#2980B9", inherit.aes = FALSE) +
  facet_wrap(~CellType, scales = "free") +
  theme_minimal() +
  labs(title    = "QP vs NNLS Deconvolution Comparison",
       subtitle = "Red dashed line = perfect agreement  |  Blue annotation = Pearson r",
       x = "QP Proportion (primary)", y = "NNLS Proportion (reference)")

ggsave(config$paths$method_comparison, p_comp, width = 10, height = 4, dpi = 300)
cat("Method comparison plot saved:", config$paths$method_comparison, "\n")
