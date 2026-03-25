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

### R packages
| Package      | Purpose                        |
|-------------|--------------------------------|
| `quadprog`  | QP-based deconvolution         |
| `ggplot2`   | Visualization                  |
| `survival`  | Survival analysis utilities    |
| `renv`      | Reproducible R environment     |
| `yaml`      | Read config.yaml               |
| `rmarkdown` | Report rendering               |
| `knitr`     | Report knitting                |

### Python scripts (`python/`)

| Script                   | Purpose                                                              |
|-------------------------|----------------------------------------------------------------------|
| `run_pipeline.py`       | Orchestrate the full pipeline: runs R scripts 01–05 in order         |
| `deconvolution_alt.py`  | Alternative deconvolution: NNLS (supervised) + NMF (unsupervised)   |

Install Python dependencies:
```bash
pip install -r python/requirements.txt
```

Run the full pipeline:
```bash
python python/run_pipeline.py
```

Run specific steps only:
```bash
python python/run_pipeline.py --steps 1 2
```

Run alternative deconvolution:
```bash
python python/deconvolution_alt.py
```

**NNLS vs NMF:**
- **NNLS** (Non-Negative Least Squares): supervised, reference-based — Python equivalent of the QP approach in `02_deconvolution.R`. Outputs `data/processed/cell_proportions_nnls.csv`.
- **NMF** (Non-negative Matrix Factorization): unsupervised, reference-free — discovers cell type components from data alone. Useful when a validated signature matrix is unavailable. Outputs `data/processed/cell_proportions_nmf.csv`.

## Reproducibility

Restore the exact R package versions used in this project with:

```r
renv::restore()
```

All thresholds and file paths are centralized in `config.yaml`. No values are hardcoded in the analysis scripts.
