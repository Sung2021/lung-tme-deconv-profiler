# =============================================================================
# 04_outlier_detection.R
# Detect TME composition outliers within each Stage group
# =============================================================================

set.seed(42)
library(ggplot2)

# --- Step 1: Load configuration ---
config    <- yaml::read_yaml("config.yaml")
threshold <- config$thresholds$outlier_percentile

# --- Step 2: Load cell proportions and clinical metadata ---
props    <- read.csv(config$paths$cell_proportions, check.names = FALSE)
# NOTE: metadata must contain SampleID and Stage columns
# metadata <- read.csv("data/raw/clinical_metadata.csv", check.names = FALSE)
# merged   <- merge(props, metadata[, c("SampleID", "Stage")], by = "SampleID")
merged <- props  # placeholder until metadata is available

cell_types <- c("Cancer", "Immune", "Stroma")

# --- Step 3: Per-Stage Z-scores and Mahalanobis distance ---
merged$Mahalanobis <- NA
merged$OutlierFlag <- FALSE

stages <- sort(unique(merged$Stage))

for (s in stages) {
  idx  <- which(merged$Stage == s)
  X    <- as.matrix(merged[idx, cell_types])

  # Z-scores per cell type within this stage
  mu   <- colMeans(X)
  sdev <- apply(X, 2, sd)

  for (ct in cell_types) {
    z_col <- paste0("Z_", ct)
    if (!(z_col %in% names(merged))) merged[[z_col]] <- NA
    merged[idx, z_col] <- (merged[idx, ct] - mu[ct]) / sdev[ct]
  }

  # Mahalanobis distance (Cancer + Immune + Stroma jointly)
  cov_mat <- cov(X)
  # Add small ridge for numerical stability
  cov_mat <- cov_mat + diag(1e-6, ncol(cov_mat))
  md <- mahalanobis(X, center = mu, cov = cov_mat)
  merged$Mahalanobis[idx] <- md

  # Flag outliers exceeding the percentile threshold
  cutoff <- quantile(md, probs = threshold / 100)
  merged$OutlierFlag[idx] <- md > cutoff
}

# --- Step 4: Label hidden high-risk patients ---
merged$RiskLabel <- ifelse(
  merged$OutlierFlag & merged$Stage == 1,
  "Hidden_High_Risk",
  "Normal"
)
cat("Outlier summary:\n")
print(table(merged$Stage, merged$OutlierFlag))
cat("\nHidden High-Risk (Stage 1 outliers):", sum(merged$RiskLabel == "Hidden_High_Risk"), "\n")

# --- Step 5: Save flagged results ---
dir.create(dirname(config$paths$outliers), recursive = TRUE, showWarnings = FALSE)
write.csv(merged, config$paths$outliers, row.names = FALSE)
cat("Outlier-flagged data saved to:", config$paths$outliers, "\n")

# --- Step 6: Scatter plot — Immune vs Stroma with outliers highlighted ---
p <- ggplot(merged, aes(x = Immune, y = Stroma, color = OutlierFlag)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  facet_wrap(~Stage, labeller = label_both) +
  theme_minimal() +
  labs(
    title = "Immune vs Stroma Proportions (Outliers Highlighted)",
    x     = "Immune Proportion",
    y     = "Stroma Proportion",
    color = "Outlier"
  )

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/outlier_scatter.png", p, width = 10, height = 8, dpi = 300)
cat("Scatter plot saved to: results/figures/outlier_scatter.png\n")
