## Cross-omics partial-correlation network and Louvain clustering template.
## Public template: no local paths, raw data, participant identifiers, or private file names.

set.seed(42)

required_pkgs <- c("igraph")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages(library(igraph))

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
MANIFEST_FILE <- get_arg("manifest", Sys.getenv("NETWORK_MANIFEST_FILE", unset = NA))
OUT_DIR <- get_arg("out", Sys.getenv("OUT_DIR", unset = "results/cross_omics_network"))
ANALYSIS_CLASSES <- get_arg("analysis-classes", Sys.getenv("ANALYSIS_CLASSES", unset = ""))
COVARIATES <- get_arg("covariates", Sys.getenv("COVARIATES", unset = ""))

NODE_FDR_CUTOFF <- as.numeric(get_arg("node-fdr", Sys.getenv("NODE_FDR_CUTOFF", unset = "0.05")))
EDGE_FDR_CUTOFF <- as.numeric(get_arg("edge-fdr", Sys.getenv("EDGE_FDR_CUTOFF", unset = "0.05")))
RHO_CUTOFF <- as.numeric(get_arg("rho-cutoff", Sys.getenv("RHO_CUTOFF", unset = "0.45")))
TOP_N_PER_OMICS <- as.integer(get_arg("top-n", Sys.getenv("TOP_N_PER_OMICS", unset = "100")))
MIN_COMPLETE_N <- as.integer(get_arg("min-complete-n", Sys.getenv("MIN_COMPLETE_N", unset = "20")))
LOUVAIN_RESOLUTION <- as.numeric(get_arg("resolution", Sys.getenv("LOUVAIN_RESOLUTION", unset = "0.7")))
MIN_COMPONENT_SIZE <- as.integer(get_arg("min-component-size", Sys.getenv("MIN_COMPONENT_SIZE", unset = "1")))
SCAN_RESOLUTIONS <- get_arg("scan-resolutions", Sys.getenv("SCAN_RESOLUTIONS", unset = "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0"))

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

analysis_classes <- split_csv(ANALYSIS_CLASSES)
covariates <- split_csv(COVARIATES)
scan_resolutions <- suppressWarnings(as.numeric(split_csv(SCAN_RESOLUTIONS)))
scan_resolutions <- scan_resolutions[is.finite(scan_resolutions)]

if (is.na(PHENOTYPE_FILE) || !file.exists(PHENOTYPE_FILE)) {
  stop("Provide a phenotype CSV with --phenotype or PHENOTYPE_FILE.")
}
if (is.na(MANIFEST_FILE) || !file.exists(MANIFEST_FILE)) {
  stop("Provide a manifest CSV with --manifest or NETWORK_MANIFEST_FILE.")
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_required_csv <- function(path, label) {
  if (is.na(path) || path == "" || !file.exists(path)) stop("Missing ", label, " file: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

make_node_id <- function(omics, feature_id) {
  paste(omics, feature_id, sep = "::")
}

read_feature_matrix <- function(path) {
  df <- read_required_csv(path, "feature matrix")
  if (ncol(df) < 2) stop("Feature matrix must contain feature IDs plus sample columns: ", path)
  feature_id <- make.unique(as.character(df[[1]]))
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- feature_id
  mat
}

read_annotation <- function(path) {
  if (is.na(path) || path == "" || !file.exists(path)) return(NULL)
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("feature_id", "display_name") %in% names(df))) {
    warning("Ignoring annotation without feature_id/display_name columns: ", path)
    return(NULL)
  }
  df[, c("feature_id", "display_name"), drop = FALSE]
}

phenotype <- read_required_csv(PHENOTYPE_FILE, "phenotype")
if (!"sample_id" %in% names(phenotype)) stop("Phenotype CSV must contain sample_id.")
if (!"class_label" %in% names(phenotype)) stop("Phenotype CSV must contain class_label.")
if (anyDuplicated(phenotype$sample_id)) stop("sample_id values must be unique.")
if (length(covariates) > 0 && !all(covariates %in% names(phenotype))) {
  stop("Missing covariates in phenotype CSV: ",
       paste(setdiff(covariates, names(phenotype)), collapse = ", "))
}

phenotype$sample_id <- as.character(phenotype$sample_id)
phenotype$class_label <- as.character(phenotype$class_label)
if (length(analysis_classes) > 0) {
  phenotype <- phenotype[phenotype$class_label %in% analysis_classes, , drop = FALSE]
}
if (nrow(phenotype) < MIN_COMPLETE_N) {
  stop("Too few phenotype rows after subsetting: ", nrow(phenotype))
}

manifest <- read_required_csv(MANIFEST_FILE, "manifest")
required_manifest_cols <- c("omics", "matrix", "differential", "feature_column")
if (!all(required_manifest_cols %in% names(manifest))) {
  stop("Manifest must contain: ", paste(required_manifest_cols, collapse = ", "))
}
if (!"annotation" %in% names(manifest)) manifest$annotation <- ""

select_candidate_nodes <- function(diff_df, feature_column, omics_name, annotation) {
  required <- c(feature_column, "logFC", "P.Value", "adj.P.Val")
  if (!all(required %in% names(diff_df))) {
    stop("Differential table for ", omics_name, " must contain: ",
         paste(required, collapse = ", "))
  }
  sig <- diff_df[diff_df$adj.P.Val < NODE_FDR_CUTOFF, required, drop = FALSE]
  sig <- sig[order(sig$adj.P.Val, sig$P.Value), , drop = FALSE]
  if (nrow(sig) > TOP_N_PER_OMICS) sig <- sig[seq_len(TOP_N_PER_OMICS), , drop = FALSE]
  names(sig)[names(sig) == feature_column] <- "feature_id"
  sig$feature_id <- as.character(sig$feature_id)
  sig$omics <- omics_name
  sig$node_id <- make_node_id(omics_name, sig$feature_id)
  sig$display_name <- sig$feature_id
  if (!is.null(annotation)) {
    mapped <- annotation$display_name[match(sig$feature_id, annotation$feature_id)]
    sig$display_name <- ifelse(!is.na(mapped) & mapped != "", mapped, sig$display_name)
  }
  sig[, c("node_id", "feature_id", "display_name", "omics", "logFC", "P.Value", "adj.P.Val"),
      drop = FALSE]
}

layers <- vector("list", nrow(manifest))
names(layers) <- manifest$omics
candidate_nodes <- list()
expression_by_layer <- list()

for (i in seq_len(nrow(manifest))) {
  omics_name <- manifest$omics[[i]]
  mat <- read_feature_matrix(manifest$matrix[[i]])
  missing_samples <- setdiff(phenotype$sample_id, colnames(mat))
  if (length(missing_samples) > 0) {
    stop(omics_name, " matrix is missing sample IDs, e.g. ",
         paste(head(missing_samples, 5), collapse = ", "))
  }
  mat <- mat[, phenotype$sample_id, drop = FALSE]
  annotation <- read_annotation(manifest$annotation[[i]])
  diff_df <- read_required_csv(manifest$differential[[i]], paste(omics_name, "differential result"))
  nodes <- select_candidate_nodes(diff_df, manifest$feature_column[[i]], omics_name, annotation)

  nodes <- nodes[!duplicated(nodes$node_id), , drop = FALSE]
  candidate_nodes[[omics_name]] <- nodes
  expression_by_layer[[omics_name]] <- mat
}

node_meta <- do.call(rbind, candidate_nodes)
rownames(node_meta) <- NULL
write.csv(node_meta, file.path(OUT_DIR, "cross_omics_candidate_nodes.csv"), row.names = FALSE)

if (nrow(node_meta) < 2) stop("Fewer than two candidate nodes were selected.")
if (length(unique(node_meta$omics)) < 2) stop("At least two omics layers are required for cross-omics edges.")

expression_by_node <- list()
for (i in seq_len(nrow(node_meta))) {
  omics_name <- node_meta$omics[[i]]
  feature_id <- node_meta$feature_id[[i]]
  expression_by_node[[node_meta$node_id[[i]]]] <- as.numeric(expression_by_layer[[omics_name]][feature_id, ])
}

test_partial_correlation <- function(node_a, node_b) {
  data <- data.frame(
    x = expression_by_node[[node_a]],
    y = expression_by_node[[node_b]],
    phenotype[, covariates, drop = FALSE],
    check.names = FALSE
  )
  data <- data[stats::complete.cases(data), , drop = FALSE]
  if (nrow(data) < MIN_COMPLETE_N) return(NULL)
  result <- tryCatch(
    {
      if (length(covariates) == 0) {
        ct <- suppressWarnings(stats::cor.test(data$x, data$y, method = "spearman", exact = FALSE))
        list(estimate = unname(ct$estimate), p.value = ct$p.value)
      } else {
        ranked <- as.data.frame(lapply(data, function(col) rank(col, ties.method = "average")))
        covar_formula <- stats::as.formula(paste("value ~", paste(covariates, collapse = " + ")))
        rx <- stats::residuals(stats::lm(covar_formula, data = transform(ranked, value = x)))
        ry <- stats::residuals(stats::lm(covar_formula, data = transform(ranked, value = y)))
        ct <- suppressWarnings(stats::cor.test(rx, ry, method = "pearson"))
        list(estimate = unname(ct$estimate), p.value = ct$p.value)
      }
    },
    error = function(e) NULL
  )
  if (is.null(result)) return(NULL)
  data.frame(
    from = node_a,
    to = node_b,
    rho = as.numeric(result$estimate),
    p_value = as.numeric(result$p.value),
    n_complete = nrow(data),
    stringsAsFactors = FALSE
  )
}

omics_pairs <- combn(unique(node_meta$omics), 2, simplify = FALSE)
edge_tests <- list()
for (pair in omics_pairs) {
  nodes_a <- node_meta$node_id[node_meta$omics == pair[[1]]]
  nodes_b <- node_meta$node_id[node_meta$omics == pair[[2]]]
  for (node_a in nodes_a) {
    for (node_b in nodes_b) {
      edge_tests[[length(edge_tests) + 1]] <- test_partial_correlation(node_a, node_b)
    }
  }
}

all_edges <- do.call(rbind, edge_tests)
if (is.null(all_edges) || nrow(all_edges) == 0) {
  stop("No valid cross-omics correlation tests were completed.")
}
all_edges$FDR <- p.adjust(all_edges$p_value, method = "BH")
all_edges$correlation_type <- ifelse(all_edges$rho > 0, "Positive", "Negative")
all_edges$weight <- abs(all_edges$rho)
write.csv(all_edges, file.path(OUT_DIR, "cross_omics_partial_correlations_all.csv"), row.names = FALSE)

edges <- all_edges[all_edges$FDR < EDGE_FDR_CUTOFF & abs(all_edges$rho) > RHO_CUTOFF, , drop = FALSE]
if (nrow(edges) == 0) {
  stop("No edges passed the selected FDR and correlation thresholds.")
}
write.csv(edges, file.path(OUT_DIR, "cross_omics_edges.csv"), row.names = FALSE)

keep_nodes <- unique(c(edges$from, edges$to))
nodes_final <- node_meta[node_meta$node_id %in% keep_nodes, , drop = FALSE]

graph <- igraph::graph_from_data_frame(
  d = edges[, c("from", "to", "weight", "rho", "FDR", "correlation_type"), drop = FALSE],
  vertices = transform(nodes_final, name = node_id),
  directed = FALSE
)

component_info <- igraph::components(graph)
component_size <- component_info$csize[component_info$membership]
names(component_size) <- names(component_info$membership)
if (MIN_COMPONENT_SIZE > 1) {
  nodes_to_drop <- names(component_size)[component_size < MIN_COMPONENT_SIZE]
  if (length(nodes_to_drop) > 0) graph <- igraph::delete_vertices(graph, nodes_to_drop)
}
if (igraph::vcount(graph) == 0 || igraph::ecount(graph) == 0) {
  stop("Network became empty after component filtering.")
}

if (length(scan_resolutions) > 0) {
  scan <- do.call(rbind, lapply(scan_resolutions, function(resolution_value) {
    set.seed(42)
    clustering <- igraph::cluster_louvain(
      graph,
      weights = igraph::E(graph)$weight,
      resolution = resolution_value
    )
    data.frame(
      resolution = resolution_value,
      n_modules = length(unique(igraph::membership(clustering))),
      module_sizes = paste(as.integer(table(igraph::membership(clustering))), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(scan, file.path(OUT_DIR, "louvain_resolution_scan.csv"), row.names = FALSE)
}

set.seed(42)
louvain <- igraph::cluster_louvain(
  graph,
  weights = igraph::E(graph)$weight,
  resolution = LOUVAIN_RESOLUTION
)

igraph::V(graph)$module <- as.integer(igraph::membership(louvain))
igraph::V(graph)$degree <- igraph::degree(graph)

node_output <- data.frame(
  node_id = igraph::V(graph)$name,
  module = igraph::V(graph)$module,
  degree = igraph::V(graph)$degree,
  stringsAsFactors = FALSE
)
node_output <- merge(node_output, node_meta, by = "node_id", all.x = TRUE)
node_output <- node_output[order(node_output$module, -node_output$degree, node_output$omics), ]

edge_output <- edges[edges$from %in% node_output$node_id & edges$to %in% node_output$node_id, , drop = FALSE]
edge_output$from_name <- node_meta$display_name[match(edge_output$from, node_meta$node_id)]
edge_output$to_name <- node_meta$display_name[match(edge_output$to, node_meta$node_id)]
edge_output <- edge_output[order(-abs(edge_output$rho), edge_output$FDR), ]

module_summary <- aggregate(
  node_id ~ module + omics,
  data = node_output,
  FUN = length
)
names(module_summary)[names(module_summary) == "node_id"] <- "n_nodes"
module_summary <- module_summary[order(module_summary$module, module_summary$omics), ]

write.csv(node_output, file.path(OUT_DIR, "cross_omics_nodes.csv"), row.names = FALSE)
write.csv(edge_output, file.path(OUT_DIR, "cross_omics_edges_with_names.csv"), row.names = FALSE)
write.csv(module_summary, file.path(OUT_DIR, "cross_omics_module_summary.csv"), row.names = FALSE)

message("Completed cross-omics network analysis.")
message("Nodes: ", nrow(node_output), "; edges: ", nrow(edge_output),
        "; Louvain modules: ", length(unique(node_output$module)))
message("Outputs written to: ", normalizePath(OUT_DIR, mustWork = FALSE))
