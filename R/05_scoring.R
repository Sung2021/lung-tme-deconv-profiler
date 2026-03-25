# =============================================================================
# 05_scoring.R
# Compute Immune Response Score and Recurrence Risk Score per patient
# =============================================================================

set.seed(42)

# --- Step 1: Load configuration ---
config <- yaml::read_yaml("config.yaml")
immune_cutoff <- config$thresholds$immune_score_cutoff
stroma_ratio  <- config$thresholds$stroma_risk_ratio

# --- Step 2: Load cell proportions and stage profiles ---
props          <- read.csv(config$paths$cell_proportions, check.names = FALSE)
stage_profiles <- read.csv(config$paths$stage_profiles, check.names = FALSE)
outlier_data   <- read.csv(config$paths$outliers, check.names = FALSE)

# --- Step 3: Compute Immune Response Score ---
# ImmuneScore = Immune / (Immune + Cancer)
props$ImmuneScore <- props$Immune / (props$Immune + props$Cancer)
props$ImmuneFlag  <- ifelse(props$ImmuneScore < immune_cutoff, "Low", "Adequate")
cat("Immune Response Score computed.\n")

# --- Step 4: Compute Recurrence Risk Score ---
# RecurrenceScore = Patient Stroma% / Stage-mean Stroma%
# Get stage-mean Stroma for each patient's Stage
stroma_means <- stage_profiles[stage_profiles$CellType == "Stroma", c("Stage", "Mean")]
names(stroma_means) <- c("Stage", "StromaStageMean")

props <- merge(props, outlier_data[, c("SampleID", "Stage", "OutlierFlag", "RiskLabel")], by = "SampleID")
props <- merge(props, stroma_means, by = "Stage")

props$RecurrenceScore <- props$Stroma / props$StromaStageMean
props$RecurrenceFlag  <- ifelse(props$RecurrenceScore > stroma_ratio, "High", "Normal")
cat("Recurrence Risk Score computed.\n")

# --- Step 5: Assemble final patient report table ---
report <- props[, c(
  "SampleID", "Stage",
  "ImmuneScore", "ImmuneFlag",
  "RecurrenceScore", "RecurrenceFlag",
  "OutlierFlag", "RiskLabel"
)]
names(report)[names(report) == "RiskLabel"] <- "RiskSummary"

dir.create(dirname(config$paths$patient_scores), recursive = TRUE, showWarnings = FALSE)
write.csv(report, config$paths$patient_scores, row.names = FALSE)
cat("Patient scores saved to:", config$paths$patient_scores, "\n")

# --- Step 6: Print single-sample example summary ---
cat("\n========== Example Patient Report Card ==========\n")
example <- report[1, ]
cat("SampleID:         ", example$SampleID, "\n")
cat("Stage:            ", example$Stage, "\n")
cat("Immune Score:     ", round(example$ImmuneScore, 3), " [", example$ImmuneFlag, "]\n")
cat("Recurrence Score: ", round(example$RecurrenceScore, 3), " [", example$RecurrenceFlag, "]\n")
cat("Outlier Flag:     ", example$OutlierFlag, "\n")
cat("Risk Summary:     ", example$RiskSummary, "\n")
cat("==================================================\n")
