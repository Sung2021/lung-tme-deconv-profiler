# TME-Resolution: Cellular Composition Index for Lung Cancer

## Background

Traditional TNM staging classifies lung cancer patients based on tumor size, lymph node involvement, and distant metastasis. However, TNM staging alone does not capture the cellular heterogeneity within the **Tumor Microenvironment (TME)**, which plays a critical role in immunotherapy response and disease recurrence. This project uses bulk RNA-seq deconvolution to estimate TME cellular composition (Cancer, Immune, Stroma fractions) and derives two clinically actionable scores: an **Immune Response Score** predicting immunotherapy responsiveness, and a **Recurrence Risk Score** identifying hidden high-risk patients within early-stage (Stage 1) lung cancer.

## Pipeline Flowchart

```
Raw counts → Normalization → QP Deconvolution → Stage Profiling → Outlier Detection → Scoring → Report
```

## How to Run

1. Restore the R environment:
   ```r
   renv::restore()
   ```
2. Execute R scripts in order:
   ```
   Rscript R/01_preprocessing.R
   Rscript R/02_deconvolution.R
   Rscript R/03_stage_profiling.R
   Rscript R/04_outlier_detection.R
   Rscript R/05_scoring.R
   ```
3. Render the report:
   ```r
   rmarkdown::render("report/analysis_report.Rmd")
   ```

## Dependencies

| Package      | Purpose                        |
|-------------|--------------------------------|
| `quadprog`  | QP-based deconvolution         |
| `ggplot2`   | Visualization                  |
| `survival`  | Survival analysis utilities    |
| `renv`      | Reproducible R environment     |
| `yaml`      | Read config.yaml               |
| `rmarkdown` | Report rendering               |
| `knitr`     | Report knitting                |

## Reproducibility

Restore the exact R package versions used in this project with:

```r
renv::restore()
```

All thresholds and file paths are centralized in `config.yaml`. No values are hardcoded in the analysis scripts.
