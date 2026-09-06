# Cross-Omics Network and Louvain Clustering Template

This folder contains a cleaned, data-free template for cross-omics network integration and Louvain community detection. The workflow builds a network from cross-omics partial correlations among differential features and exports node, edge, and module tables for downstream visualization.

The template does not contain local paths, raw data, participant identifiers, or private project file names.

## Workflow

1. Select candidate nodes from each omics layer using differential-analysis FDR.
2. Extract feature values from selected samples, such as two disease groups or an analysis subset.
3. Test pairwise cross-omics partial Spearman correlations while adjusting for user-specified covariates.
4. Retain edges passing both edge-level FDR and absolute correlation thresholds.
5. Build an undirected weighted network with `igraph`.
6. Optionally remove small connected components before final module detection.
7. Run Louvain community detection and export node, edge, and module tables.

## Input Manifest

Create a CSV manifest with one row per omics layer:

```csv
omics,matrix,differential,feature_column,annotation
Metabolomics,data/metabolomics_matrix.csv,data/metab_target_vs_reference.csv,metabolite,
Lipidomics,data/lipidomics_matrix.csv,data/lipid_target_vs_reference.csv,lipid,
Proteomics,data/proteomics_matrix.csv,data/prot_target_vs_reference.csv,protein,data/protein_annotation.csv
```

Required formats:

- `matrix`: first column is `feature_id`; remaining columns are samples.
- `differential`: contains the feature ID column named by `feature_column`, plus `logFC`, `P.Value`, and `adj.P.Val`.
- `annotation`, optional: contains `feature_id` and `display_name`.

The phenotype CSV must contain:

- `sample_id`: sample identifiers matching matrix columns.
- `class_label`: group labels.
- Any covariates named in `--covariates`, such as `age`, `sex`, or another study-specific adjustment set.

## Example

```bash
Rscript build_cross_omics_louvain_network_template.R \
  --phenotype data/phenotype.csv \
  --manifest data/network_manifest.csv \
  --analysis-classes HM,UM \
  --covariates age,sex,iop_mean \
  --node-fdr 0.05 \
  --edge-fdr 0.05 \
  --rho-cutoff 0.45 \
  --top-n 100 \
  --resolution 0.7 \
  --out results/cross_omics_network
```

## Outputs

- `cross_omics_candidate_nodes.csv`: candidate differential features.
- `cross_omics_partial_correlations_all.csv`: all tested cross-omics partial correlations.
- `cross_omics_edges.csv`: retained network edges.
- `cross_omics_nodes.csv`: retained nodes with Louvain module and degree.
- `cross_omics_module_summary.csv`: module size summaries.
- `louvain_resolution_scan.csv`: optional module counts across tested resolutions.

Raw participant-level data and local analysis outputs are intentionally excluded from this repository.
