# =============================================================================
# 01_preprocessing.R
# Normalize raw RNA-seq count matrix to log2 CPM
# =============================================================================

set.seed(42)

# --- Step 1: Load configuration ---
config <- yaml::read_yaml("config.yaml")

# --- Step 2: Load raw count matrix ---
counts <- read.csv(config$paths$raw_counts, row.names = 1, check.names = FALSE)
cat("Raw count matrix loaded:", nrow(counts), "genes x", ncol(counts), "samples\n")

# --- Step 3: Filter low-expression genes ---
keep <- rowSums(counts) > 10
counts_filtered <- counts[keep, ]
cat("Genes after filtering (rowSums > 10):", nrow(counts_filtered), "\n")

# --- Step 4: Log2 CPM normalization ---
lib_sizes <- colSums(counts_filtered)
cpm <- sweep(counts_filtered, 2, lib_sizes, FUN = "/") * 1e6
log2cpm <- log2(cpm + 1)
cat("Log2 CPM normalization complete.\n")

# --- Step 5: Save normalized expression matrix ---
dir.create(dirname(config$paths$normalized_expr), recursive = TRUE, showWarnings = FALSE)
write.csv(log2cpm, config$paths$normalized_expr, row.names = TRUE)
cat("Normalized expression saved to:", config$paths$normalized_expr, "\n")
