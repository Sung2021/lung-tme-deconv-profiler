# Architecture

## Overview

**lung-tme-deconv-profiler** is a modular, reproducible pipeline for decomposing bulk RNA-seq to estimate tumor microenvironment (TME) cellular composition and derive clinically actionable scoring metrics for lung cancer patients.

## Project Layout

```
lung-tme-deconv-profiler/
├── R/                 # Ordered R stages and the 02 wrapper shim
├── python/            # Orchestrator and optional alternate deconvolution
├── data/
│   ├── raw/           # Input counts and sample metadata
│   ├── reference/     # Signature matrices
│   └── processed/     # Normalized matrices and cell proportions
├── results/
│   ├── tables/        # Stage profiles, outliers, patient scores
│   └── figures/       # QC and comparison plots
├── report/            # Renderable analysis report
└── config.yaml        # Single source of truth for paths and thresholds
```

### Input Contract

The pipeline assumes three inputs before Stage 1 starts:

- `data/raw/counts.csv`: raw gene-by-sample count matrix with genes in rows and samples in columns
- `data/reference/sig_matrix.csv`: reference signature matrix with the same gene ID type and normalization space as the expression matrix
- `data/raw/*metadata*.csv`: sample metadata containing at minimum `SampleID` and `Stage`

The current R implementation still contains placeholders in [R/03_stage_profiling.R](R/03_stage_profiling.R) and [R/04_outlier_detection.R](R/04_outlier_detection.R), so the metadata source must be made explicit before those stages can be run end-to-end.

### Validation Gates

Each stage should be treated as passing only if its local checks succeed:

- Stage 1: gene retention after filtering is sufficient and TPM normalization completes
- Stage 2a: signature matrix has non-zero variance and acceptable condition number
- Stage 2b: QP and NNLS agree closely across cell types
- Stage 3: each Stage group has enough samples to estimate mean and SD reliably
- Stage 4: stage-wise covariance is invertible or ridge-stabilized before Mahalanobis scoring
- Stage 5: patient scores are produced without missing joins across proportions, stage profiles, and outlier flags

### Pipeline Flow Diagram

```
Raw Counts (counts.csv)
        ↓
[01] Preprocessing: Normalize to log₂ TPM
        ↓ (expr_normalized.csv)
[02a] Signature Matrix: Select robustly expressed cell-type markers
        ↓ (sig_matrix_validated.csv)
[02b] Deconvolution: QP + NNLS solve ||S·w - y||² → cell proportions
        ↓ (proportions_qp.csv, proportions_nnls.csv)
[03] Stage Profiling: Compute mean cellular composition per TNM stage
        ↓ (stage_profiles.csv)
[04] Outlier Detection: Flag samples with extreme Mahalanobis distance
        ↓ (outliers_flagged.csv)
[05] Scoring: Compute Immune Response & Recurrence Risk scores per patient
        ↓ (patient_scores.csv, analysis_report.html)
Final Clinical Report
```

---

## Detailed Stage Description

### **Stage 1: Preprocessing (01_preprocessing.R)**

**Purpose:** Normalize raw RNA-seq counts to log₂ TPM (Transcripts Per Million).

**Key Steps:**
1. Load raw count matrix from `config$paths$raw_counts`
2. **Dual filter:**
   - Keep genes expressed (count > 0) in ≥70% of samples (prevalence cutoff)
   - Keep genes with `rowSum ≥ 14` (minimum aggregate signal)
3. Fetch gene lengths from Ensembl BioMart using `biomaRt` package
   - Uses longest transcript per gene as representative length
   - Configured via `config$biomart$dataset` and `config$biomart$gene_id_type`
4. Compute TPM: `normalized[i, j] = log₂(1 + (count[i, j] / gene_length[i]) × 1e6 / library_size[j])`
   - TPM is sample-comparable (each sample sums to ~1e6 before log transform)
   - Corrects for gene length bias (unlike CPM)

**Outputs:**
- `data/processed/expr_normalized.csv` — subset of genes × all samples, TPM space

**Why TPM?**
- Comparable across samples (library-size normalized)
- Corrects gene-length bias (essential for deconvolution)
- Preferred by signature matrix builders (e.g., immunedeconv, xCell)

---

### **Stage 2a: Signature Matrix Validation (02a_signature_matrix.R)**

**Purpose:** Identify and validate a set of robustly expressed, cell-type-specific marker genes.

**Key Steps:**
1. Load reference signature (a curated set of cell-type markers)
   - Expected format: genes × cell types (Cancer, Immune, Stroma)
2. Intersect with normalized expression; filter to available genes
3. Compute per-cell-type statistics (mean, variance, specificity)
4. **Robustness check:**
   - Require markers with low variance and high mean expression across reference samples
   - Flag collinear or redundant markers
5. Output validated signature matrix: only robust markers retained

**Outputs:**
- `data/reference/sig_matrix_validated.csv` — filtered signature (typically 10–50 genes × 3 cell types)

**Why Validation?**
- Prevents numerical instability in the deconvolution solver (Stage 2b)
- Improves biological interpretability (strong, specific markers per cell type)

---

### **Stage 2b: Deconvolution (02b_deconvolution.R)**

**Purpose:** Estimate cell-type proportions (Cancer, Immune, Stroma) for each sample via constrained linear algebra.

**Mathematical Formulation:**
```
Given:
  S ∈ ℝ^(m × k)  — validated signature matrix (m genes, k cell types)
  y ∈ ℝ^m        — sample expression vector

Solve:
  min ||S·w - y||²₂  (least-squares fit of signature basis to sample)
  w∈ℝᵏ

Subject to:
  w_i ≥ 0        for all i  (non-negative cell proportions)
  Σ w_i = 1               (proportions sum to 1 — biological hard constraint)
```

**Two Methods Implemented:**

| Method | Solver | Constraints | Role |
|--------|--------|-------------|------|
| **QP** *(Primary)* | Quadratic Programming (R `quadprog`) | `w_i ≥ 0` AND `Σw_i = 1` enforced simultaneously | Canonical output; enforces sum-to-one as hard constraint |
| **NNLS** *(Reference)* | Non-Negative Least Squares | `w_i ≥ 0` only; sum-to-one applied post-hoc | Sanity check; if QP & NNLS agree (r > 0.95 per cell type), results are numerically robust |

**Why QP is Primary:**
1. Sum-to-one is a biological hard constraint, not a post-hoc correction
2. Enforced during optimization → avoids normalization errors
3. Mathematical basis of commercial CIBERSORT (uses SVR + similar constraints)

**Outputs:**
- `data/processed/proportions_qp.csv` — QP results
- `data/processed/proportions_nnls.csv` — NNLS results (comparison)
- `data/processed/cell_proportions.csv` — canonical QP output (used downstream)
- `results/figures/cell_composition_barplot.png` — stacked barplot: samples × cell types
- `results/figures/method_comparison.png` — QP vs. NNLS scatter per cell type

**Entry Point Note:**
- `R/02_deconvolution.R` is a wrapper kept for backward compatibility
- The canonical execution path is `02a_signature_matrix.R -> 02b_deconvolution.R`

---

### **Stage 3: Stage Profiling (03_stage_profiling.R)**

**Purpose:** Establish reference TME composition profiles stratified by TNM stage.

**Key Steps:**
1. Load cell proportions (from Stage 2b)
2. Join sample metadata containing `SampleID` and `Stage`
3. For each stage and cell type, compute:
   - Mean proportion
   - Std. deviation
   - Median, 95th percentile (for outlier detection threshold setting)
4. Build reference profile table: `Stage × CellType → Mean, SD, P95`

**Outputs:**
- `results/tables/stage_profiles.csv` 
  - Columns: `Stage, CellType, Mean, SD, P95`
  - Used in Stage 4 (outlier detection) and Stage 5 (recurrence scoring)

**Clinical Insight:**
Early-stage (Stage 1) tumors with high stroma burden (exceeding Stage 1 mean) may indicate hidden disease progression risk despite TNM classification.

---

### **Stage 4: Outlier Detection (04_outlier_detection.R)**

**Purpose:** Flag samples with anomalous cellular composition (using Mahalanobis distance).

**Key Steps:**
1. Load cell proportions (3-dimensional vector: Cancer, Immune, Stroma)
2. Join Stage metadata before grouping by stage
3. For each stage group, compute covariance matrix of proportions
4. Compute Mahalanobis distance: `D²_M = (x - μ)ᵀ Σ⁻¹ (x - μ)`
   - Accounts for correlations between cell types
   - More robust than Euclidean distance
5. Set outlier threshold: `D²_M > percentile_p` (default: p = 95th percentile, configurable via `config$thresholds$outlier_percentile`)
6. Flag samples with extreme values; assign risk labels

**Outputs:**
- `results/tables/outliers_flagged.csv`
  - Columns: `SampleID, Stage, Mahalanobis_Distance, IsOutlier, RiskLabel`
  - `RiskLabel`: e.g., "Extreme_Immune" (very high immune), "Extreme_Stroma" (very high stroma), "Normal"

**Clinical Use:**
Outlier samples often indicate:
- Biological heterogeneity (mixed histology, necrosis, immune infiltration)
- Hidden disease progression (e.g., high stroma in Stage 1 specimens)

---

### **Stage 5: Scoring (05_scoring.R)**

**Purpose:** Compute two clinically actionable summary scores per patient.

#### **Score 1: Immune Response Score**

```
ImmuneScore = Immune / (Immune + Cancer)  ∈ [0, 1]
```

- Quantifies the balance of immune cells to malignant cells
- **Biological interpretation:** High score → favorable immune infiltration → better immunotherapy response potential
- **Clinical cutoff (default):** `ImmuneScore ≥ 0.20` (configurable: `config$thresholds$immune_score_cutoff`)
  - `ImmuneScore < 0.20` → flagged as "Low" (poor immune infiltration)
  - `ImmuneScore ≥ 0.20` → flagged as "Adequate"

**Predictive Hypothesis:**
Samples with adequate immune score demonstrate higher checkpoint inhibitor response rates.

---

#### **Score 2: Recurrence Risk Score**

```
RecurrenceScore = Patient_Stroma% / Stage_Mean_Stroma%
```

- Compares a sample's stroma burden to the mean stroma burden of its TNM stage
- **Biological interpretation:** High stroma in early-stage tumors may indicate hidden high-risk features (desmoplasia, fibroblast activation, immunosuppression)
- **Clinical cutoff (default):** `RecurrenceScore > 1.50` (configurable: `config$thresholds$stroma_risk_ratio`)
  - `RecurrenceScore > 1.50` → flagged as "High" (elevated recurrence risk)
  - `RecurrenceScore ≤ 1.50` → flagged as "Normal"

**Predictive Hypothesis:**
Stage 1 patients with high stroma burden (RecurrenceScore > 1.50) show elevated *actual* recurrence risk despite early-stage classification, identifying patients who may benefit from adjuvant therapy.

---

**Outputs:**
- `results/tables/patient_scores.csv`
  - Columns: `SampleID, Stage, ImmuneScore, ImmuneFlag, RecurrenceScore, RecurrenceFlag, OutlierFlag, RiskSummary`
  - Each row = one patient's clinical summary
- Generated into `report/analysis_report.Rmd` as dynamic table + figures

---

## Configuration

All hardcoded values are centralized in **`config.yaml`**:

```yaml
thresholds:
  immune_score_cutoff: 0.20          # Immune / (Immune + Cancer) threshold
  stroma_risk_ratio: 1.50            # RecurrenceScore threshold
  outlier_percentile: 95             # Mahalanobis distance cutoff

biomart:
  dataset: hsapiens_gene_ensembl     # Ensembl dataset (change for non-human)
  gene_id_type: hgnc_symbol          # Gene ID format in count matrix

paths:
  raw_counts: data/raw/counts.csv
  signature_matrix: data/reference/sig_matrix.csv
  sig_matrix_validated: data/reference/sig_matrix_validated.csv
  normalized_expr: data/processed/expr_normalized.csv
  cell_proportions: data/processed/cell_proportions.csv
  proportions_qp: data/processed/proportions_qp.csv
  proportions_nnls: data/processed/proportions_nnls.csv
  stage_profiles: results/tables/stage_profiles.csv
  outliers: results/tables/outliers_flagged.csv
  patient_scores: results/tables/patient_scores.csv
  method_comparison: results/figures/method_comparison.png
```

**Usage:**
- Each R script runs `yaml::read_yaml("config.yaml")` at start; no hardcoding
- Python orchestrator (`python/run_pipeline.py`) also respects `config.yaml` for reproducibility
- **To customize:** Edit `config.yaml` once; all stages adapt

## Execution Order

The canonical pipeline is:

1. `R/01_preprocessing.R`
2. `R/02a_signature_matrix.R`
3. `R/02b_deconvolution.R`
4. `R/03_stage_profiling.R`
5. `R/04_outlier_detection.R`
6. `R/05_scoring.R`

`R/02_deconvolution.R` remains available only as a compatibility wrapper for older calls.

---

## Data Flow & Dependency Graph

```
config.yaml
    ↓
data/raw/counts.csv ──→ [01_preprocessing] ──→ expr_normalized.csv
                                                    ↓
data/reference/sig_matrix.csv ──→ [02a_signature_matrix] ──→ sig_matrix_validated.csv
                                        ↓
                                    [02b_deconvolution] ──→ proportions_qp.csv
                                        ↓                    proportions_nnls.csv
                                        ↓                    cell_composition_barplot.png
                                        ↓                    method_comparison.png
                                        ↓
                                  [03_stage_profiling] ──→ stage_profiles.csv
                                        ↓
                                  [04_outlier_detection] ──→ outliers_flagged.csv
                                        ↓
                                    [05_scoring] ──→ patient_scores.csv
                                        ↓
                                  [Report Rendering] ──→ analysis_report.html
```

---

## Reproducibility & Environment

### R Environment
```r
renv::restore()  # Restores exact package versions from renv.lock
```

### Python Environment
```bash
pip install -r python/requirements.txt
```

### Running the Pipeline

**Option 1: Sequential R execution**
```bash
Rscript R/01_preprocessing.R
Rscript R/02a_signature_matrix.R
Rscript R/02b_deconvolution.R
Rscript R/03_stage_profiling.R
Rscript R/04_outlier_detection.R
Rscript R/05_scoring.R
rmarkdown::render("report/analysis_report.Rmd")
```

**Option 2: Python orchestrator (recommended)**
```bash
# Run all steps
python python/run_pipeline.py

# Run specific steps only
python python/run_pipeline.py --steps 1 2 5
```

**Option 3: Alternative Deconvolution (NNLS + NMF)**
```bash
python python/deconvolution_alt.py
```

---

## Key Design Decisions

### 1. **Why Log₂ TPM?**
- **Not RPKM:** RPKM is not comparable across samples; TPM is normalized to constant total (1e6 per sample)
- **Not CPM:** CPM ignores gene length; TPM corrects for it
- **Not raw counts:** Deconvolution requires normalized, comparable values

### 2. **Why QP + NNLS Comparison?**
- **QP:** Enforces sum-to-one as hard constraint → mathematically principled
- **NNLS:** Post-hoc sum-to-one normalization → sanity check for numerical stability
- **Agreement (r > 0.95)** = evidence of robustness; **disagreement** = re-examine signature

### 3. **Why Mahalanobis Distance?**
- **Not Euclidean:** Euclidean ignores covariance structure; Mahalanobis accounts for correlations between cell types
- **Multivariate outlier detection:** More robust when proportions are correlated (e.g., high Immune ↔ low Stroma)

### 4. **Why Two Scores?**
- **Immune Response Score:** Predicts immunotherapy response (T-cell infiltration)
- **Recurrence Risk Score:** Identifies hidden high-risk early-stage patients (desmoplasia, immune exclusion)
- **Complementary:** One captures infiltration, the other captures suppressive stromal context

### 5. **Why Config-Driven?**
- No hardcoding → reproducible, versioned, audit-able
- Easy sensitivity analysis (adjust thresholds, re-run)
- Multi-dataset adaptability (change Ensembl dataset, gene ID type, etc.)

---

## Output Summary

| Stage | Output File | Format | Purpose |
|-------|-------------|--------|---------|
| 1 | `expr_normalized.csv` | CSV (genes × samples) | Normalized expression matrix (input to deconvolution) |
| 2a | `sig_matrix_validated.csv` | CSV (genes × cell types) | Validated signature for deconvolution solver |
| 2b | `proportions_qp.csv` | CSV (samples × cell types) | QP cell proportions (canonical) |
| 2b | `proportions_nnls.csv` | CSV (samples × cell types) | NNLS cell proportions (sanity check) |
| 2b | `cell_composition_barplot.png` | PNG | Stacked barplot of cell composition |
| 2b | `method_comparison.png` | PNG | Scatter: QP vs. NNLS per cell type |
| 3 | `stage_profiles.csv` | CSV | Reference profiles (Stage × CellType) |
| 4 | `outliers_flagged.csv` | CSV | Outlier samples + Mahalanobis distance |
| 5 | `patient_scores.csv` | CSV | Clinical summary per patient |
| After 5 | `analysis_report.html` | HTML (interactive) | Full reproducible report (Rmarkdown output) |

---

## Extensions & Future Directions

1. **Spatial deconvolution:** Adapt for spatial transcriptomics (Visium, MERFISH)
2. **Alternative deconvolution methods:** Immunedeconv, Seurat, CellDART for comparison
3. **Prognostic validation:** Prospective survival analysis on independent cohort
4. **Multi-omics integration:** Include TCR/BCR diversity, surface protein expression
5. **Docker containerization:** Ensure reproducibility across HPC/cloud environments
