# =============================================================================
# 01_preprocessing.R
# Normalize raw RNA-seq count matrix to log2 TPM
# Gene lengths are fetched automatically from Ensembl via biomaRt.
# TPM is preferred over CPM for deconvolution: it corrects for gene length
# bias and produces sample-comparable proportions (each sample sums to 1e6).
# =============================================================================

set.seed(42)

# --- Step 1: Load configuration ---
config <- yaml::read_yaml("config.yaml")

# --- Step 2: Load raw count matrix ---
counts <- read.csv(config$paths$raw_counts, row.names = 1, check.names = FALSE)
cat("Raw count matrix loaded:", nrow(counts), "genes x", ncol(counts), "samples\n")

# --- Step 3: Filter low-expression genes (dual cutoff) ---
# Condition 1: count > 0 in at least 70% of samples (expression prevalence)
prev_pass  <- rowMeans(counts > 0) >= 0.70
# Condition 2: total row sum >= 1*2*7 = 14 (minimum aggregate signal)
sum_pass   <- rowSums(counts) >= 1 * 2 * 7
keep       <- prev_pass & sum_pass
counts_filtered <- counts[keep, ]
cat(sprintf(
  "Genes after filtering (prevalence>=70%% AND rowSums>=%d): %d (removed %d)\n",
  1L * 2L * 7L, nrow(counts_filtered), nrow(counts) - nrow(counts_filtered)
))

# --- Step 4: Fetch gene lengths from Ensembl via biomaRt ---
cat("Connecting to Ensembl BioMart (", config$biomart$dataset, ")...\n")
mart <- biomaRt::useMart("ensembl", dataset = config$biomart$dataset)
gene_info <- biomaRt::getBM(
  attributes = c(config$biomart$gene_id_type, "transcript_length"),
  filters    = config$biomart$gene_id_type,
  values     = rownames(counts_filtered),
  mart       = mart
)
cat("Gene length records retrieved:", nrow(gene_info), "\n")

# Use the longest transcript per gene as its representative length
gene_lengths <- tapply(
  gene_info[["transcript_length"]],
  gene_info[[config$biomart$gene_id_type]],
  max
)

# Retain only genes for which a length was retrieved
common_genes    <- intersect(rownames(counts_filtered), names(gene_lengths))
counts_tpm      <- as.matrix(counts_filtered[common_genes, ])
lengths_kb      <- gene_lengths[common_genes] / 1000
cat(sprintf("Genes retained after length matching: %d (dropped %d)\n",
            length(common_genes), nrow(counts_filtered) - length(common_genes)))

# --- Step 5: TPM normalization ---
# Step 5a: RPK — reads per kilobase of transcript
rpk     <- sweep(counts_tpm, 1, lengths_kb, FUN = "/")
# Step 5b: per-sample scaling factor (sum of RPK / 1e6)
scaling <- colSums(rpk) / 1e6
# Step 5c: TPM
tpm     <- sweep(rpk, 2, scaling, FUN = "/")
log2tpm <- log2(tpm + 1)
cat("Log2 TPM normalization complete. Sample TPM column sums (should all equal 1e6):\n")
print(round(colSums(tpm)[1:min(5, ncol(tpm))], 0))  # sanity-check first 5 samples

# --- Step 6: Save normalized expression matrix ---
dir.create(dirname(config$paths$normalized_expr), recursive = TRUE, showWarnings = FALSE)
write.csv(log2tpm, config$paths$normalized_expr, row.names = TRUE)
cat("Normalized expression saved to:", config$paths$normalized_expr, "\n")
