# Ultra-Long Axial Myopia Multi-Omics

This repository contains reproducible analysis code for a plasma multi-omics study of ultra-long axial myopia. The project integrates lipidomic, metabolomic, and proteomic profiles to characterize systemic molecular alterations across normal control, high myopia, and ultra-long axial myopia groups.

The study aims to describe multi-layer molecular remodeling associated with extreme axial elongation, identify severity-related molecular patterns, and explore cross-omics signals linked to lipid metabolism, mitochondrial and energy-related pathways, neuro-glial injury, and inflammatory responses. The analysis workflow includes differential omics analysis, pathway and functional enrichment, multi-omics integration, feature selection, diagnostic model development, and manuscript figure generation.

This codebase is intended to support transparent and reproducible reporting of the computational workflow. Raw participant-level data are not included in this repository and should be accessed only through the approved project data environment.

## Analysis Code

- `analysis/lasso_multiclass/`: cleaned template scripts for single-omics and multi-omics multiclass LASSO modeling, including repeated cross-validation, fold-specific feature prefiltering, model evaluation, feature stability summaries, and clinical-anchor plotting.
