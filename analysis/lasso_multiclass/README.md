# Multiclass LASSO Template

This folder contains cleaned, data-free template scripts for multiclass LASSO modeling in plasma multi-omics data. The scripts were prepared for public code sharing and do not contain local paths, raw data, participant identifiers, or project-specific private file names.

## Input Format

The main script expects one phenotype table and one or more omics matrices.

Phenotype table:

- CSV file with one row per sample.
- Required columns: `sample_id`, `class_label`.
- `class_label` should contain the analysis groups in the intended order, for example `NC`, `HM`, and `UM`.

Omics matrices:

- CSV files with features in rows and samples in columns.
- The first column must contain feature identifiers.
- Remaining column names must match `sample_id` in the phenotype table.

Optional protein annotation:

- CSV file with columns `feature_id` and `display_name`.
- Used only to improve feature labels in exported coefficient tables.

## Main Analysis

Example:

```bash
Rscript run_multiclass_lasso_template.R \
  --phenotype data/phenotype.csv \
  --metabolomics data/metabolomics_matrix.csv \
  --lipidomics data/lipidomics_matrix.csv \
  --proteomics data/proteomics_matrix.csv \
  --protein-annotation data/protein_annotation.csv \
  --classes NC,HM,UM \
  --out results/lasso_multiclass
```

The workflow fits separate single-omics models and an integrated multi-omics model. It performs repeated stratified outer cross-validation, conducts univariate limma prefiltering inside each training fold, applies fold-specific imputation and scaling, fits multinomial LASSO models with `glmnet`, and exports performance metrics, out-of-fold probabilities, feature stability, and final display-model coefficients.

## Clinical Anchor Plot

Use `plot_signature_score_vs_clinical_grade_template.R` to relate a selected class probability, such as the UM signature score, to an external ordered clinical grade.

Example:

```bash
Rscript plot_signature_score_vs_clinical_grade_template.R \
  --predictions results/lasso_multiclass/Multiomics/Multiomics_sample_level_predictions.csv \
  --clinical data/clinical_grade.csv \
  --target-class UM \
  --grade-column clinical_grade \
  --out results/lasso_multiclass/clinical_anchor
```

## Notes

The scripts are templates for reproducible analysis and manuscript reporting. Raw participant-level data and local project outputs are intentionally excluded.
