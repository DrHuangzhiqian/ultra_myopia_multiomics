# Molecular Pattern Reclassification Template

This folder contains a cleaned, data-free template for reclassifying multi-omics features into mutually exclusive molecular patterns across three ordered groups, such as control, intermediate disease, and advanced disease.

The template was derived from the Figure 4 molecular-pattern workflow and prepared for public code sharing. It does not contain local paths, raw data, participant identifiers, or private project file names.

## Pattern Logic

The script prioritizes group-specific patterns before residual severity-gradient patterns:

- `Target-specific`: features significant in target vs control, but not significant in intermediate vs control.
- `Intermediate-specific`: features significant in intermediate vs control, with disappearance in target vs control or reversal in target vs intermediate.
- `Gradient`: remaining features associated with the ordered phenotype severity score by Spearman correlation.

This priority order enforces mutually exclusive feature assignment.

## Input Manifest

Create a CSV manifest with one row per omics layer:

```csv
omics,matrix,diff_target_control,diff_intermediate_control,diff_target_intermediate,correlation,annotation
Lipidomics,data/lipidomics_matrix.csv,data/lipid_target_control.csv,data/lipid_intermediate_control.csv,data/lipid_target_intermediate.csv,data/lipid_severity_correlation.csv,
Metabolomics,data/metabolomics_matrix.csv,data/metab_target_control.csv,data/metab_intermediate_control.csv,data/metab_target_intermediate.csv,data/metab_severity_correlation.csv,
Proteomics,data/proteomics_matrix.csv,data/prot_target_control.csv,data/prot_intermediate_control.csv,data/prot_target_intermediate.csv,data/prot_severity_correlation.csv,data/protein_annotation.csv
```

Required table formats:

- `matrix`: first column is `feature_id`; remaining columns are samples.
- differential tables: `feature_id`, `logFC`, `P.Value`, `adj.P.Val`.
- correlation table: `feature_id`, `rho`, `p_value`.
- annotation table, optional: `feature_id`, `display_name`.

The phenotype CSV must contain:

- `sample_id`: sample identifiers matching matrix columns.
- `class_label`: group labels, such as `NC`, `HM`, and `UM`.
- `severity_score`: numeric ordered score used for column ordering and optional checks.

## Example

```bash
Rscript reclassify_molecular_patterns_template.R \
  --phenotype data/phenotype.csv \
  --manifest data/pattern_manifest.csv \
  --classes NC,HM,UM \
  --control-class NC \
  --intermediate-class HM \
  --target-class UM \
  --out results/molecular_patterns
```

## Outputs

- `pattern_candidate_features.csv`: all candidate features assigned to each pattern.
- `pattern_counts.csv`: feature counts by omics layer and pattern.
- `pattern_heatmap_features.csv`: selected features shown in heatmaps.
- `pattern_heatmap_grid.pdf/png`: combined multi-panel figure.
- `target_specific_heatmap.pdf/png`, `intermediate_specific_heatmap.pdf/png`, `gradient_heatmap.pdf/png`: separate heatmaps.
- `pattern_counts_barplot.pdf/png`: count summary.

Raw participant-level data and local analysis outputs are intentionally excluded from this repository.
