## Reclassify multi-omics features into mutually exclusive molecular patterns.
## Public template: no local paths, raw data, or participant identifiers.

set.seed(42)

required_pkgs <- c("dplyr", "ggplot2", "patchwork", "pheatmap", "grid")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(grid)
})

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
MANIFEST_FILE <- get_arg("manifest", Sys.getenv("PATTERN_MANIFEST_FILE", unset = NA))
OUT_DIR <- get_arg("out", Sys.getenv("OUT_DIR", unset = "results/molecular_patterns"))
CLASS_LEVELS <- strsplit(get_arg("classes", Sys.getenv("CLASS_LEVELS", unset = "NC,HM,UM")),
                         ",", fixed = TRUE)[[1]]
CLASS_LEVELS <- trimws(CLASS_LEVELS)
CONTROL_CLASS <- get_arg("control-class", CLASS_LEVELS[[1]])
INTERMEDIATE_CLASS <- get_arg("intermediate-class", CLASS_LEVELS[[2]])
TARGET_CLASS <- get_arg("target-class", CLASS_LEVELS[[length(CLASS_LEVELS)]])

FDR_CUTOFF <- as.numeric(get_arg("fdr-cutoff", Sys.getenv("FDR_CUTOFF", unset = "0.05")))
RHO_CUTOFF <- as.numeric(get_arg("rho-cutoff", Sys.getenv("RHO_CUTOFF", unset = "0.3")))
COR_P_CUTOFF <- as.numeric(get_arg("cor-p-cutoff", Sys.getenv("COR_P_CUTOFF", unset = "0.05")))
N_UP <- as.integer(get_arg("n-up", Sys.getenv("N_UP", unset = "8")))
N_DOWN <- as.integer(get_arg("n-down", Sys.getenv("N_DOWN", unset = "8")))

if (is.na(PHENOTYPE_FILE) || !file.exists(PHENOTYPE_FILE)) {
  stop("Provide a phenotype CSV with --phenotype or PHENOTYPE_FILE.")
}
if (is.na(MANIFEST_FILE) || !file.exists(MANIFEST_FILE)) {
  stop("Provide a manifest CSV with --manifest or PATTERN_MANIFEST_FILE.")
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_required_csv <- function(path, label) {
  if (is.na(path) || path == "" || !file.exists(path)) {
    stop("Missing ", label, " file: ", path)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

read_matrix <- function(path) {
  df <- read_required_csv(path, "feature matrix")
  if (ncol(df) < 2) stop("Feature matrix must include feature IDs plus sample columns: ", path)
  feature_id <- make.unique(as.character(df[[1]]))
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- feature_id
  mat
}

standardize_diff_table <- function(path) {
  df <- read_required_csv(path, "differential result")
  required <- c("feature_id", "logFC", "P.Value", "adj.P.Val")
  if (!all(required %in% names(df))) {
    stop("Differential result must contain: ", paste(required, collapse = ", "), ". File: ", path)
  }
  df[, required, drop = FALSE]
}

standardize_correlation_table <- function(path) {
  df <- read_required_csv(path, "correlation result")
  required <- c("feature_id", "rho", "p_value")
  if (!all(required %in% names(df))) {
    stop("Correlation result must contain: ", paste(required, collapse = ", "), ". File: ", path)
  }
  df[, required, drop = FALSE]
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
if (!all(c("sample_id", "class_label", "severity_score") %in% names(phenotype))) {
  stop("Phenotype CSV must contain sample_id, class_label, and severity_score.")
}
phenotype$sample_id <- as.character(phenotype$sample_id)
phenotype$class_label <- factor(as.character(phenotype$class_label), levels = CLASS_LEVELS)
phenotype$severity_score <- as.numeric(phenotype$severity_score)
if (anyNA(phenotype$class_label)) stop("Some phenotype class labels are absent from --classes.")
if (anyDuplicated(phenotype$sample_id)) stop("sample_id values must be unique.")

manifest <- read_required_csv(MANIFEST_FILE, "manifest")
required_manifest_cols <- c(
  "omics", "matrix", "diff_target_control", "diff_intermediate_control",
  "diff_target_intermediate", "correlation"
)
if (!all(required_manifest_cols %in% names(manifest))) {
  stop("Manifest must contain: ", paste(required_manifest_cols, collapse = ", "))
}
if (!"annotation" %in% names(manifest)) manifest$annotation <- ""

omics_layers <- lapply(seq_len(nrow(manifest)), function(i) {
  omics_name <- manifest$omics[[i]]
  mat <- read_matrix(manifest$matrix[[i]])
  missing_samples <- setdiff(phenotype$sample_id, colnames(mat))
  if (length(missing_samples) > 0) {
    stop(omics_name, " matrix is missing sample IDs, e.g. ",
         paste(head(missing_samples, 5), collapse = ", "))
  }
  mat <- mat[, phenotype$sample_id, drop = FALSE]
  list(
    omics = omics_name,
    matrix = mat,
    target_control = standardize_diff_table(manifest$diff_target_control[[i]]),
    intermediate_control = standardize_diff_table(manifest$diff_intermediate_control[[i]]),
    target_intermediate = standardize_diff_table(manifest$diff_target_intermediate[[i]]),
    correlation = standardize_correlation_table(manifest$correlation[[i]]),
    annotation = read_annotation(manifest$annotation[[i]])
  )
})
names(omics_layers) <- manifest$omics
OMICS_ORDER <- names(omics_layers)

select_top_by_direction <- function(df, n_up = N_UP, n_down = N_DOWN, order_col = "P.Value") {
  if (nrow(df) == 0) return(list(up = character(0), down = character(0)))
  df <- df[order(df[[order_col]], decreasing = FALSE), , drop = FALSE]
  up <- df[df$logFC > 0, , drop = FALSE]
  down <- df[df$logFC < 0, , drop = FALSE]
  list(up = head(as.character(up$feature_id), n_up),
       down = head(as.character(down$feature_id), n_down))
}

feature_set <- function(pattern_obj) {
  unique(as.character(pattern_obj$candidates$feature_id))
}

pick_target_specific <- function(target_control, intermediate_control) {
  target_sig <- target_control[target_control$adj.P.Val < FDR_CUTOFF, , drop = FALSE]
  intermediate_fdr <- intermediate_control[, c("feature_id", "adj.P.Val"), drop = FALSE]
  colnames(intermediate_fdr) <- c("feature_id", "fdr_intermediate_control")
  merged <- merge(target_sig, intermediate_fdr, by = "feature_id", all.x = TRUE)
  merged <- merged[is.na(merged$fdr_intermediate_control) |
                     merged$fdr_intermediate_control >= FDR_CUTOFF, , drop = FALSE]
  merged <- merged[order(merged$P.Value), , drop = FALSE]
  selected <- select_top_by_direction(merged)
  candidates <- data.frame(
    feature_id = as.character(merged$feature_id),
    pattern = "Target-specific",
    logFC = merged$logFC,
    P.Value = merged$P.Value,
    adj.P.Val = merged$adj.P.Val,
    stringsAsFactors = FALSE
  )
  list(up = selected$up, down = selected$down, candidates = candidates)
}

pick_intermediate_specific <- function(intermediate_control, target_control, target_intermediate,
                                       exclude_features = character(0)) {
  intermediate_sig <- intermediate_control[intermediate_control$adj.P.Val < FDR_CUTOFF, , drop = FALSE]
  intermediate_sig <- intermediate_sig[!intermediate_sig$feature_id %in% exclude_features, , drop = FALSE]

  target_status <- target_control[, c("feature_id", "adj.P.Val", "logFC"), drop = FALSE]
  colnames(target_status) <- c("feature_id", "fdr_target_control", "logFC_target_control")
  merged <- merge(intermediate_sig, target_status, by = "feature_id", all.x = TRUE)

  target_vs_intermediate <- target_intermediate[, c("feature_id", "logFC"), drop = FALSE]
  colnames(target_vs_intermediate) <- c("feature_id", "logFC_target_intermediate")
  merged <- merge(merged, target_vs_intermediate, by = "feature_id", all.x = TRUE)

  disappeared <- is.na(merged$fdr_target_control) | merged$fdr_target_control >= FDR_CUTOFF
  reversed <- !is.na(merged$logFC_target_intermediate) &
    sign(merged$logFC) != sign(merged$logFC_target_intermediate)
  merged <- merged[disappeared | reversed, , drop = FALSE]
  merged <- merged[order(merged$P.Value), , drop = FALSE]

  selected <- select_top_by_direction(merged)
  candidates <- data.frame(
    feature_id = as.character(merged$feature_id),
    pattern = "Intermediate-specific",
    logFC = merged$logFC,
    P.Value = merged$P.Value,
    adj.P.Val = merged$adj.P.Val,
    stringsAsFactors = FALSE
  )
  list(up = selected$up, down = selected$down, candidates = candidates)
}

pick_gradient <- function(correlation, target_control, mat, exclude_features = character(0)) {
  base <- correlation[!correlation$feature_id %in% exclude_features, , drop = FALSE]
  base <- base[abs(base$rho) > RHO_CUTOFF & base$p_value < COR_P_CUTOFF, , drop = FALSE]

  diff_subset <- target_control[, c("feature_id", "logFC", "P.Value", "adj.P.Val"), drop = FALSE]
  merged <- merge(base, diff_subset, by = "feature_id", all.x = TRUE)
  merged <- merged[merged$feature_id %in% rownames(mat), , drop = FALSE]

  up <- merged[merged$rho > 0, , drop = FALSE]
  down <- merged[merged$rho < 0, , drop = FALSE]
  up <- up[order(-up$rho), , drop = FALSE]
  down <- down[order(down$rho), , drop = FALSE]

  candidates <- data.frame(
    feature_id = as.character(merged$feature_id),
    pattern = "Gradient",
    rho = merged$rho,
    p_value = merged$p_value,
    logFC = merged$logFC,
    P.Value = merged$P.Value,
    adj.P.Val = merged$adj.P.Val,
    stringsAsFactors = FALSE
  )
  list(up = head(as.character(up$feature_id), N_UP),
       down = head(as.character(down$feature_id), N_DOWN),
       candidates = candidates)
}

assign_patterns <- function(layer) {
  target <- pick_target_specific(layer$target_control, layer$intermediate_control)
  assigned <- feature_set(target)
  intermediate <- pick_intermediate_specific(
    layer$intermediate_control,
    layer$target_control,
    layer$target_intermediate,
    exclude_features = assigned
  )
  assigned <- unique(c(assigned, feature_set(intermediate)))
  gradient <- pick_gradient(
    layer$correlation,
    layer$target_control,
    layer$matrix,
    exclude_features = assigned
  )
  list(Target = target, Intermediate = intermediate, Gradient = gradient)
}

patterns <- lapply(omics_layers, assign_patterns)

display_feature <- function(feature_id, layer) {
  annotation <- layer$annotation
  if (is.null(annotation)) return(as.character(feature_id))
  label <- annotation$display_name[match(feature_id, annotation$feature_id)]
  ifelse(is.na(label) | label == "", as.character(feature_id), as.character(label))
}

build_feature_table <- function() {
  rows <- list()
  pattern_map <- c(
    Target = "Target-specific",
    Intermediate = "Intermediate-specific",
    Gradient = "Gradient"
  )
  for (omics_name in OMICS_ORDER) {
    for (pattern_key in names(pattern_map)) {
      candidates <- patterns[[omics_name]][[pattern_key]]$candidates
      if (nrow(candidates) == 0) next
      candidates$omics <- omics_name
      candidates$display_name <- vapply(
        candidates$feature_id,
        display_feature,
        character(1),
        layer = omics_layers[[omics_name]]
      )
      candidates$selected_for_heatmap <- candidates$feature_id %in%
        c(patterns[[omics_name]][[pattern_key]]$up,
          patterns[[omics_name]][[pattern_key]]$down)
      candidates$pattern <- pattern_map[[pattern_key]]
      rows[[paste(omics_name, pattern_key, sep = "_")]] <- candidates
    }
  }
  dplyr::bind_rows(rows) %>%
    dplyr::select(omics, pattern, feature_id, display_name, selected_for_heatmap, dplyr::everything())
}

feature_table <- build_feature_table()
write.csv(feature_table, file.path(OUT_DIR, "pattern_candidate_features.csv"), row.names = FALSE)

count_table <- feature_table %>%
  dplyr::count(omics, pattern, name = "n_features") %>%
  dplyr::mutate(
    omics = factor(omics, levels = OMICS_ORDER),
    pattern = factor(pattern, levels = c("Target-specific", "Intermediate-specific", "Gradient"))
  ) %>%
  dplyr::arrange(pattern, omics)
write.csv(count_table, file.path(OUT_DIR, "pattern_counts.csv"), row.names = FALSE)

selected_table <- feature_table %>%
  dplyr::filter(selected_for_heatmap) %>%
  dplyr::arrange(
    factor(pattern, levels = c("Target-specific", "Intermediate-specific", "Gradient")),
    factor(omics, levels = OMICS_ORDER)
  )
write.csv(selected_table, file.path(OUT_DIR, "pattern_heatmap_features.csv"), row.names = FALSE)

class_colors <- setNames(c("#2166AC", "#F4A900", "#C0392B")[seq_along(CLASS_LEVELS)], CLASS_LEVELS)
omics_colors <- setNames(c("#E69F00", "#009E73", "#7B2D8B", "#5B84B1", "#999999")[seq_along(OMICS_ORDER)], OMICS_ORDER)
pattern_colors <- c(
  "Target-specific" = "#C0392B",
  "Intermediate-specific" = "#F4A900",
  "Gradient" = "#5B84B1"
)

theme_pattern <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.tag = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.line = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4)
  )

make_heatmap <- function(pattern_key, title_text, tag_text, max_name_len = 30) {
  direction_order <- c("Up", "Down")
  row_info <- list()
  for (direction in direction_order) {
    for (omics_name in OMICS_ORDER) {
      pattern_obj <- patterns[[omics_name]][[pattern_key]]
      features <- if (direction == "Up") pattern_obj$up else pattern_obj$down
      features <- features[features %in% rownames(omics_layers[[omics_name]]$matrix)]
      if (length(features) == 0) next
      row_info[[paste(omics_name, direction, sep = "_")]] <- data.frame(
        feature_id = features,
        omics = omics_name,
        direction = direction,
        stringsAsFactors = FALSE
      )
    }
  }
  row_info <- dplyr::bind_rows(row_info)
  if (nrow(row_info) == 0) {
    return(grid::textGrob(paste(tag_text, "No selected features")))
  }

  expression_matrix <- do.call(rbind, lapply(seq_len(nrow(row_info)), function(i) {
    omics_layers[[row_info$omics[[i]]]]$matrix[row_info$feature_id[[i]], , drop = FALSE]
  }))
  expression_matrix <- t(scale(t(expression_matrix)))
  expression_matrix[!is.finite(expression_matrix)] <- 0

  labels <- mapply(
    function(feature_id, omics_name) display_feature(feature_id, omics_layers[[omics_name]]),
    row_info$feature_id,
    row_info$omics,
    USE.NAMES = FALSE
  )
  labels <- ifelse(nchar(labels) > max_name_len, paste0(substr(labels, 1, max_name_len), "..."), labels)
  rownames(expression_matrix) <- make.unique(labels, sep = "_")

  column_order <- order(phenotype$class_label, phenotype$severity_score, phenotype$sample_id)
  expression_matrix <- expression_matrix[, column_order, drop = FALSE]
  group_run <- rle(as.character(phenotype$class_label[column_order]))
  gaps_col <- cumsum(group_run$lengths)
  gaps_col <- gaps_col[-length(gaps_col)]

  annotation_col <- data.frame(
    Class = factor(phenotype$class_label[column_order], levels = CLASS_LEVELS),
    row.names = colnames(expression_matrix)
  )
  annotation_row <- data.frame(
    Omics = factor(row_info$omics, levels = OMICS_ORDER),
    Direction = factor(row_info$direction, levels = direction_order),
    row.names = rownames(expression_matrix)
  )
  direction_run <- rle(as.character(annotation_row$Direction))
  gaps_row <- cumsum(direction_run$lengths)
  gaps_row <- gaps_row[-length(gaps_row)]

  heatmap_obj <- pheatmap::pheatmap(
    expression_matrix,
    color = colorRampPalette(c("#2166AC", "#F7F7F7", "#C0392B"))(100),
    breaks = seq(-2.5, 2.5, length.out = 101),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = list(
      Class = class_colors,
      Omics = omics_colors,
      Direction = c(Up = "#C0392B", Down = "#2166AC")
    ),
    gaps_col = gaps_col,
    gaps_row = gaps_row,
    fontsize_row = 8,
    fontsize = 10,
    border_color = NA,
    main = paste0(tag_text, "   ", title_text),
    silent = TRUE
  )
  heatmap_obj$gtable
}

n_for_pattern <- function(pattern_label) {
  sum(count_table$n_features[count_table$pattern == pattern_label], na.rm = TRUE)
}

grob_target <- make_heatmap("Target", sprintf("Target-specific (n=%d)", n_for_pattern("Target-specific")), "A")
grob_intermediate <- make_heatmap("Intermediate", sprintf("Intermediate-specific (n=%d)", n_for_pattern("Intermediate-specific")), "B")
grob_gradient <- make_heatmap("Gradient", sprintf("Gradient (n=%d)", n_for_pattern("Gradient")), "C")

count_plot <- ggplot(count_table, aes(x = omics, y = n_features, fill = pattern)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(aes(label = n_features), position = position_dodge(width = 0.72),
            vjust = -0.45, size = 3.6, fontface = "bold") +
  scale_x_discrete(limits = OMICS_ORDER) +
  scale_fill_manual(values = pattern_colors, name = "Pattern") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Number of features", tag = "D") +
  theme_pattern +
  theme(plot.tag.position = c(0, 1), legend.position = "right",
        axis.text.x = element_text(angle = 20, hjust = 1))

combined_plot <- (patchwork::wrap_elements(grob_target) | patchwork::wrap_elements(grob_intermediate)) /
  (patchwork::wrap_elements(grob_gradient) | count_plot) +
  patchwork::plot_layout(heights = c(1, 1), widths = c(1.2, 1))

save_grob <- function(grob, filename_base, width = 9, height = 10) {
  pdf(file.path(OUT_DIR, paste0(filename_base, ".pdf")), width = width, height = height, useDingbats = FALSE)
  grid::grid.newpage()
  grid::grid.draw(grob)
  dev.off()
  png(file.path(OUT_DIR, paste0(filename_base, ".png")), width = width, height = height, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(grob)
  dev.off()
}

pdf(file.path(OUT_DIR, "pattern_heatmap_grid.pdf"), width = 20, height = 18, useDingbats = FALSE)
print(combined_plot)
dev.off()
png(file.path(OUT_DIR, "pattern_heatmap_grid.png"), width = 20, height = 18, units = "in", res = 300)
print(combined_plot)
dev.off()

save_grob(grob_target, "target_specific_heatmap")
save_grob(grob_intermediate, "intermediate_specific_heatmap")
save_grob(grob_gradient, "gradient_heatmap")
ggsave(file.path(OUT_DIR, "pattern_counts_barplot.pdf"), count_plot, width = 7, height = 6)
ggsave(file.path(OUT_DIR, "pattern_counts_barplot.png"), count_plot, width = 7, height = 6, dpi = 300)

message("Completed molecular pattern reclassification. Outputs written to: ",
        normalizePath(OUT_DIR, mustWork = FALSE))
