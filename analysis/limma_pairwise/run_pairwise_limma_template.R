## Pairwise limma differential analysis template for metabolomics, lipidomics, and proteomics.
## Public template: no local paths, raw data, participant identifiers, plotting code, or prediction-model code.

set.seed(42)

required_pkgs <- c("limma")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages(library(limma))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--")) stop("Unexpected argument: ", args[[i]])
    key <- sub("^--", "", args[[i]])
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

args <- parse_args()
get_arg <- function(name, default = NULL) {
  value <- args[[name]]
  if (is.null(value) || identical(value, "")) default else value
}

PHENOTYPE_FILE <- get_arg("phenotype", Sys.getenv("PHENOTYPE_FILE", unset = NA))
METABOLOMICS_FILE <- get_arg("metabolomics", Sys.getenv("METABOLOMICS_FILE", unset = NA))
LIPIDOMICS_FILE <- get_arg("lipidomics", Sys.getenv("LIPIDOMICS_FILE", unset = NA))
PROTEOMICS_FILE <- get_arg("proteomics", Sys.getenv("PROTEOMICS_FILE", unset = NA))
OUT_DIR <- get_arg("out", Sys.getenv("OUT_DIR", unset = "results/limma_pairwise"))

CLASS_LEVELS <- strsplit(get_arg("classes", Sys.getenv("CLASS_LEVELS", unset = "NC,HM,UM")),
                         ",", fixed = TRUE)[[1]]
CLASS_LEVELS <- trimws(CLASS_LEVELS)

FDR_CUTOFF <- as.numeric(get_arg("fdr-cutoff", Sys.getenv("FDR_CUTOFF", unset = "0.05")))
LOGFC_CUTOFF <- as.numeric(get_arg("logfc-cutoff", Sys.getenv("LOGFC_CUTOFF", unset = "0.58")))

default_comparisons <- function(class_levels) {
  if (length(class_levels) == 3) {
    return(c(
      paste(class_levels[[3]], class_levels[[1]], sep = "_vs_"),
      paste(class_levels[[3]], class_levels[[2]], sep = "_vs_"),
      paste(class_levels[[2]], class_levels[[1]], sep = "_vs_")
    ))
  }
  comparisons <- character(0)
  for (i in seq_along(class_levels)[-1]) {
    for (j in seq_len(i - 1)) {
      comparisons <- c(comparisons, paste(class_levels[[i]], class_levels[[j]], sep = "_vs_"))
    }
  }
  comparisons
}

COMPARISONS <- get_arg("comparisons", paste(default_comparisons(CLASS_LEVELS), collapse = ","))
COMPARISONS <- trimws(strsplit(COMPARISONS, ",", fixed = TRUE)[[1]])

if (is.na(PHENOTYPE_FILE) || !file.exists(PHENOTYPE_FILE)) {
  stop("Provide a phenotype CSV with --phenotype or PHENOTYPE_FILE.")
}
if (length(CLASS_LEVELS) < 2) stop("At least two class levels are required.")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_feature_matrix <- function(path) {
  if (is.na(path) || path == "" || !file.exists(path)) return(NULL)
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) < 2) stop("Feature matrix must contain feature IDs plus sample columns: ", path)
  feature_id <- make.unique(as.character(df[[1]]))
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- feature_id
  mat
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

parse_comparison <- function(comparison) {
  parts <- strsplit(comparison, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2 || any(!nzchar(parts))) {
    stop("Comparison must use the format ClassA_vs_ClassB: ", comparison)
  }
  if (!all(parts %in% CLASS_LEVELS)) {
    stop("Comparison contains labels absent from --classes: ", comparison)
  }
  list(target = parts[[1]], reference = parts[[2]])
}

make_contrast_matrix <- function(comparisons, design_colnames) {
  contrast <- matrix(0, nrow = length(design_colnames), ncol = length(comparisons),
                     dimnames = list(design_colnames, comparisons))
  for (comparison in comparisons) {
    spec <- parse_comparison(comparison)
    target <- make.names(spec$target)
    reference <- make.names(spec$reference)
    contrast[target, comparison] <- 1
    contrast[reference, comparison] <- -1
  }
  contrast
}

annotate_direction <- function(result, logfc_cutoff = LOGFC_CUTOFF, fdr_cutoff = FDR_CUTOFF) {
  result$Direction <- "No change"
  result$Direction[result$logFC > logfc_cutoff & result$adj.P.Val < fdr_cutoff] <- "Up"
  result$Direction[result$logFC < -logfc_cutoff & result$adj.P.Val < fdr_cutoff] <- "Down"
  result
}

phenotype <- read.csv(PHENOTYPE_FILE, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample_id", "class_label") %in% names(phenotype))) {
  stop("Phenotype CSV must contain columns: sample_id, class_label.")
}
phenotype$sample_id <- as.character(phenotype$sample_id)
phenotype$class_label <- factor(as.character(phenotype$class_label), levels = CLASS_LEVELS)
if (anyNA(phenotype$class_label)) {
  stop("Some class_label values are absent from --classes.")
}
if (anyDuplicated(phenotype$sample_id)) stop("sample_id values must be unique.")

omics_mats <- list(
  Metabolomics = read_feature_matrix(METABOLOMICS_FILE),
  Lipidomics = read_feature_matrix(LIPIDOMICS_FILE),
  Proteomics = read_feature_matrix(PROTEOMICS_FILE)
)
omics_mats <- omics_mats[!vapply(omics_mats, is.null, logical(1))]
if (length(omics_mats) == 0) {
  stop("Provide at least one feature matrix with --metabolomics, --lipidomics, or --proteomics.")
}

design <- model.matrix(~ 0 + phenotype$class_label)
colnames(design) <- make.names(CLASS_LEVELS)
contrast <- make_contrast_matrix(COMPARISONS, colnames(design))

run_limma_for_layer <- function(layer_name, mat) {
  missing_samples <- setdiff(phenotype$sample_id, colnames(mat))
  if (length(missing_samples) > 0) {
    stop(layer_name, " matrix is missing sample IDs, e.g. ",
         paste(head(missing_samples, 5), collapse = ", "))
  }
  mat <- mat[, phenotype$sample_id, drop = FALSE]

  layer_dir <- file.path(OUT_DIR, safe_name(layer_name))
  dir.create(layer_dir, recursive = TRUE, showWarnings = FALSE)

  fit <- limma::lmFit(mat, design)
  fit2 <- limma::contrasts.fit(fit, contrast)
  fit2 <- limma::eBayes(fit2)

  summary_rows <- list()
  for (comparison in COMPARISONS) {
    result <- limma::topTable(
      fit2,
      coef = comparison,
      number = Inf,
      sort.by = "P",
      adjust.method = "BH"
    )
    result$feature_id <- rownames(result)
    result$omics <- layer_name
    result$comparison <- comparison
    result <- annotate_direction(result)
    result <- result[, c("omics", "comparison", "feature_id",
                         setdiff(names(result), c("omics", "comparison", "feature_id"))),
                     drop = FALSE]

    full_file <- file.path(layer_dir, paste0(safe_name(layer_name), "_", safe_name(comparison), "_full.csv"))
    sig_file <- file.path(layer_dir, paste0(safe_name(layer_name), "_", safe_name(comparison), "_FDRsig.csv"))
    write.csv(result, full_file, row.names = FALSE)
    write.csv(result[result$Direction != "No change", , drop = FALSE], sig_file, row.names = FALSE)

    summary_rows[[comparison]] <- data.frame(
      omics = layer_name,
      comparison = comparison,
      n_features_total = nrow(result),
      n_up = sum(result$Direction == "Up", na.rm = TRUE),
      n_down = sum(result$Direction == "Down", na.rm = TRUE),
      n_fdr_significant = sum(result$Direction != "No change", na.rm = TRUE),
      logfc_cutoff = LOGFC_CUTOFF,
      fdr_cutoff = FDR_CUTOFF,
      stringsAsFactors = FALSE
    )
  }

  layer_summary <- do.call(rbind, summary_rows)
  write.csv(layer_summary, file.path(layer_dir, paste0(safe_name(layer_name), "_limma_summary_counts.csv")),
            row.names = FALSE)
  layer_summary
}

all_summaries <- do.call(rbind, Map(run_limma_for_layer, names(omics_mats), omics_mats))
write.csv(all_summaries, file.path(OUT_DIR, "all_omics_limma_summary_counts.csv"), row.names = FALSE)

message("Completed pairwise limma differential analysis. Outputs written to: ",
        normalizePath(OUT_DIR, mustWork = FALSE))
