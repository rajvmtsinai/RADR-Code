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

---

## Software and R Packages

R (recent version) with the following packages: survival, stats, gridExtra, broom, tidyr, stringr, dplyr, reshape2, ggplot2, ggpubr, survminer, openxlsx

If you are building the dataset directly from VCF files, you may also need: `vcfR`.

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

---

## Outputs

* **Figures:** Forest plots print during the run. To save high-resolution PDFs, uncomment the `ggsave(...)` lines in each notebook and specify a relative path within this folder (e.g., `out/plots/...`).

---

## Reproducibility Checklist
* [ ] R and all listed packages installed.
* [ ] Input files available in this folder (or subfolders) with correct column names, including `Sample/sample_id`, `GENDER`, `AGE_OF_ENTRY`, `PC1..PC10`, per-phenotype binaries, `Genotype`, and `CHROM/POS/REF/ALT`.
* [ ] Paths set relative to this folder (avoid hard-coded absolute paths).
* [ ] Minimum carrier count thresholds (if any) are set appropriately for your sample size.

---
