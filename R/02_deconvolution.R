# =============================================================================
# 02_deconvolution.R
# QP-based deconvolution to estimate cell type proportions per sample
# =============================================================================

set.seed(42)
library(quadprog)
library(ggplot2)

# --- Step 1: Load configuration ---
config <- yaml::read_yaml("config.yaml")

# --- Step 2: Load normalized expression and signature matrix ---
expr <- read.csv(config$paths$normalized_expr, row.names = 1, check.names = FALSE)
sig  <- read.csv(config$paths$signature_matrix, row.names = 1, check.names = FALSE)
cat("Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("Signature matrix:", nrow(sig), "genes x", ncol(sig), "cell types\n")

# --- Step 3: Align genes between expression and signature ---
common_genes <- intersect(rownames(expr), rownames(sig))
cat("Common genes used for deconvolution:", length(common_genes), "\n")
expr_aligned <- as.matrix(expr[common_genes, ])
sig_aligned  <- as.matrix(sig[common_genes, ])

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
