# =============================================================================
# 03_stage_profiling.R
# Compute cell type statistics per Stage and test for differences
# =============================================================================

set.seed(42)
library(ggplot2)

# --- Step 1: Load configuration ---
config <- yaml::read_yaml("config.yaml")

# --- Step 2: Load cell proportions and clinical metadata ---
props    <- read.csv(config$paths$cell_proportions, check.names = FALSE)
metadata <- read.csv(config$paths$raw_counts, check.names = FALSE)  # placeholder: replace with actual metadata path
# NOTE: metadata must contain columns SampleID and Stage (values 1–4)
# For now, assume a metadata file is loaded separately and merged
# metadata <- read.csv("data/raw/clinical_metadata.csv", check.names = FALSE)

# --- Step 3: Merge proportions with Stage info ---
# Ensure metadata has SampleID and Stage columns
merged <- merge(props, metadata[, c("SampleID", "Stage")], by = "SampleID")
cell_types <- c("Cancer", "Immune", "Stroma")

# --- Step 4: Compute mean and SD per Stage per cell type ---
stage_summary <- do.call(rbind, lapply(cell_types, function(ct) {
  do.call(rbind, lapply(sort(unique(merged$Stage)), function(s) {
    vals <- merged[merged$Stage == s, ct]
    data.frame(
      CellType = ct,
      Stage    = s,
      Mean     = mean(vals, na.rm = TRUE),
      SD       = sd(vals, na.rm = TRUE),
      N        = length(vals)
    )
  }))
}))
cat("Stage summary:\n")
print(stage_summary)

# --- Step 5: One-way ANOVA per cell type across stages ---
anova_results <- data.frame(CellType = cell_types, P_value = NA)
for (i in seq_along(cell_types)) {
  ct <- cell_types[i]
  fit <- aov(as.formula(paste(ct, "~ factor(Stage)")), data = merged)
  anova_results$P_value[i] <- summary(fit)[[1]][["Pr(>F)"]][1]
}
cat("\nANOVA results:\n")
print(anova_results)

# --- Step 6: Append ANOVA p-values to summary and save ---
stage_summary <- merge(stage_summary, anova_results, by = "CellType")
dir.create(dirname(config$paths$stage_profiles), recursive = TRUE, showWarnings = FALSE)
write.csv(stage_summary, config$paths$stage_profiles, row.names = FALSE)
cat("Stage profiles saved to:", config$paths$stage_profiles, "\n")

# --- Step 7: Boxplot of cell proportions by Stage ---
plot_df <- tidyr::pivot_longer(
  merged,
  cols      = all_of(cell_types),
  names_to  = "CellType",
  values_to = "Proportion"
)

p <- ggplot(plot_df, aes(x = factor(Stage), y = Proportion, fill = CellType)) +
  geom_boxplot() +
  facet_wrap(~CellType, scales = "free_y") +
  theme_minimal() +
  labs(title = "Cell Type Proportions by Stage", x = "Stage", y = "Proportion")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/stage_boxplot.png", p, width = 10, height = 6, dpi = 300)
cat("Boxplot saved to: results/figures/stage_boxplot.png\n")
