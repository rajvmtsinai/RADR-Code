# ADRD Variant Regression Pipeline

This README describes how to use the project files in the respository. It explains what each file does, what inputs are expected, how to run the analyses, and where outputs are produced.

---

## Project Contents

- **`Regression File Creation.Rmd`**  
  Builds a regression-ready table by merging VCF genotypes with variant annotations and cohort covariates/phenotypes. Generates a single binary outcome (e.g., `dementia_present`) for downstream models.

- **`Multiple Regression.Rmd`**  
  For each variant, fits a multivariable logistic regression (e.g., `dementia_present ~ carrier + GENDER + AGE_OF_ENTRY + PC1..PC10`), extracts odds ratios (OR) with 95% CIs, applies FDR correction, and draws a log-scale forest plot.  
  A knitted HTML copy may be present as: **`Multiple Regression.nb.html`**

- **`Cox Regression.Rmd`**  
  For each variant, fits a Cox proportional hazards model on time-to-onset (e.g., `Surv(AGE_OF_ENTRY, dementia_present) ~ carrier + GENDER + PCs`), extracts hazard ratios (HR) with 95% CIs, and draws a log-scale forest plot.  
  A knitted HTML copy may be present as: **`Cox Regression.nb.html`**

- **`PC_Plots.Rmd`**  
  Creates ancestry principal component (PC) scatter plots for the full cohort and for carriers of each P/LP/RF variant. It merges PC coordinates with genetic ancestry labels and sample IDs, identifies variant carriers from VCF genotypes, and generates per-variant PC1 vs PC2 visualizations with consistent ancestry color mapping.

- **`APOE_cox_logistic_regression.Rmd`**
  Runs per-variant APOE-adjusted Cox and logistic regression models from a regression-ready CSV plus an APOE sample table. It standardizes APOE labels, derives E2/E4 dosage, filters variants by minimum carrier count, and writes APOE-adjusted effect-size result tables (including optional interaction tests).

---

## Software and R Packages

R (recent version) with the following packages: survival, stats, gridExtra, broom, tidyr, stringr, dplyr, reshape2, ggplot2, ggpubr, survminer, openxlsx, readr, purrr

If you are building the dataset directly from VCF files, you may also need: `vcfR`.

For principal-component ancestry plots, you may also need: `viridisLite`, `readxl`, `MatchIt`, and `RColorBrewer`.

---

## Expected Inputs

1. **VCF** with P/LP/RF variants and sample genotypes coded as `0/0`, `0/1`, `1/1`.

2. **Variant annotation table** (e.g., RADR export) with `CHROM/POS/REF/ALT`, `Gene`, and protein change to create readable variant labels.

3. **Phenotype table** keyed by sample ID with ADRD-related binary columns (e.g., `EOAD`, `LOAD`, `DLB`, `VD`, etc.). These are typically collapsed into a single binary outcome (`dementia_present`).  
   *Essentially this is a phenotype file with `sample_id` and associated `phenotype` as the two columns. One sample can have multiple phenotypes, which will correspond to multiple rows in the phenotype table.*

4. **Covariate table** keyed by `Sample` (or `sample_id`) including `GENDER/SEX`, `AGE_OF_ENTRY`, and principal components (`PC1` → `PC10`).  
   **Note:** `AGE_OF_ENTRY` is the age of the visit where they were diagnosed with the ADRD phenotype.

5. **(Optional)** **APOE status** table keyed by sample ID.

---

## What Each Notebook Does

### A) `Regression File Creation.Rmd`
- Reads VCF and annotation; merges on `CHROM/POS/REF/ALT` with RADR.
- Creates a human-readable variant label (`Gene:Protein_Change`).
- Reshapes genotypes to long format; defines `carrier = (Genotype != "0/0")`. *Carriers are not differentiated between heterozygous and homozygous.*
- Merges covariates and phenotypes; collapses phenotypes to `dementia_present`.  
  If they have an ADRD phenotype → `dementia_present = 1`; otherwise `0`.
- **Output:** an in-memory data frame ready for regression (optionally saved as CSV if you uncomment the write line).

### B) `Multiple Regression.Rmd` (Odds Ratios)
Model:
```r
glm(dementia_present ~ carrier + GENDER + AGE_OF_ENTRY + PC1 + PC2 + ... + PC10,
    family = binomial)
```
Extracts OR, 95% CI, p-value; applies FDR correction.

Produces a log-scale forest plot with a reference line at OR = 1.

Optional: uncomment ggsave(...) lines to write high-DPI PDF(s) to a chosen folder.

### C) Cox Regression.Rmd (Hazard Ratios)

You can add or remove APOE status to any of the models.
```r
coxph(Surv(AGE_OF_ENTRY, dementia_present) ~ carrier + GENDER + PCs)
```
Extracts HR, 95% CI, p-value; creates a log-scale forest plot with a reference at HR = 1.

Optional: uncomment ggsave(...) lines to write high-DPI PDF(s) to a chosen folder.

### D) PC_Plots.Rmd (Ancestry Visualization)
- Loads principal component coordinates and genetic ancestry calls per sample.
- Merges participant IDs with cohort metadata to align sample identifiers.
- Reads variant calls from VCF and links them to RADR annotations.
- Identifies carrier samples (`0/1` or `1/1`) for each variant.
- Creates:
  - a cohort-level PC1 vs PC2 ancestry plot, and
  - per-variant PC1 vs PC2 carrier ancestry plots.
- Uses fixed ancestry categories (`AFR`, `EAS`, `EUR`, `HIS`, `OTHER`, `SAS`) and a consistent custom palette.

Optional: uncomment `ggsave(...)` lines to export PNG files for each per-variant plot and the cohort plot.

### E) APOE_cox_logistic_regression.Rmd (APOE-Adjusted Variant Models)
- Reads a per-variant regression dataset (`cox_regression_data.csv`) and APOE table (`apoe_status.csv`) using configurable column names.
- Joins APOE onto the analysis dataset and harmonizes APOE labels to `E2/E3/E4` combinations (e.g., `E3-E3`, `E2-E4`, `E4-E4`).
- Creates APOE dosage covariates (`e2_dosage`, `e4_dosage`) and binary `carrier` from genotype (`0/0` vs carrier).
- Filters variants by minimum carrier count (`MIN_CARRIERS`) before modeling.
- Fits per-variant APOE-adjusted Cox models and saves:
  - `out/apoe_adjusted_models/cox_variant_effects_apoe_adjusted.csv`
  - `out/apoe_adjusted_models/sanity_apoe_status_counts.csv`
  - `out/apoe_adjusted_models/sanity_carrier_by_apoe.csv`
- Optionally tests carrier-by-APOE interaction (default: `carrier * e4_dosage`; optional factor interaction).
- Fits APOE-adjusted logistic models and saves:
  - `out/apoe_adjusted_models/logit_variant_effects_apoe_adjusted.csv`
  - `out/apoe_adjusted_models/significant_variants_effect_sizes_apoe_adjusted.csv`


## How to Run (Plain-Folder Setup)

### 1) Build the dataset

* Open `Regression File Creation.Rmd` in RStudio.
* Update any file paths to point to your local project folder.
* If you prefer a saved file, uncomment the `write.csv(...)` line (if present) and choose a relative path within this folder.
* The final output will be a long-type dataframe saved as CSV that can be used for the regression analysis.

### 2) Run logistic models (ORs)

* Open `Multiple Regression.Rmd`.
* Ensure the regression-ready dataset from step 1 is available in your R session or load the saved CSV.
* Run to produce results and the forest plot.

### 3) Run Cox models (HRs)

* Open `Cox Regression.Rmd`.
* Ensure the merged dataset and (optional) APOE status are available (or load them from saved CSVs).
* Run to produce results and the forest plot.

### 4) Run ancestry PC plots

* Open `PC_Plots.Rmd`.
* Update file paths for:
  - the PC file,
  - the ancestry labels file,
  - cohort/sample ID mapping,
  - VCF input,
  - RADR annotation,
  - phenotype table.
* Run all chunks to generate per-variant carrier ancestry plots and the full-cohort ancestry PC plot.

### 5) Run APOE-adjusted Cox + logistic models

* Open `APOE_cox_logistic_regression.Rmd`.
* Update the configuration block at the top for:
  - input file paths (`COX_DATA_PATH`, `APOE_DATA_PATH`),
  - sample ID columns (`COX_SAMPLE_COL`, `APOE_SAMPLE_COL`),
  - covariate names (`SEX_COL`, `PC_COLS`),
  - outcome/time columns (`EVENT_COL`, `TIME_COL`), and
  - minimum carrier threshold (`MIN_CARRIERS`).
* Confirm `TIME_COL` is age-at-event for cases and age-at-censoring for controls.
* Run all chunks to generate APOE-adjusted model outputs under `out/apoe_adjusted_models/`.

---

## Outputs

* **Figures:** Forest plots print during the run. To save high-resolution PDFs, uncomment the `ggsave(...)` lines in each notebook and specify a relative path within this folder (e.g., `out/plots/...`).
* **Ancestry plots:** `PC_Plots.Rmd` prints PC scatter plots to the plotting device. To save files, uncomment the `ggsave(...)` lines and set an output directory.

---

## Reproducibility Checklist
* [ ] R and all listed packages installed.
* [ ] Input files available in this folder (or subfolders) with correct column names, including `Sample/sample_id`, `GENDER`, `AGE_OF_ENTRY`, `PC1..PC10`, per-phenotype binaries, `Genotype`, and `CHROM/POS/REF/ALT`.
* [ ] Paths set relative to this folder (avoid hard-coded absolute paths).
* [ ] Minimum carrier count thresholds (if any) are set appropriately for your sample size.

---
