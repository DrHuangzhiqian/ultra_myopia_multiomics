## Multiclass LASSO template for single-omics and multi-omics classification.
## Public template: no local paths, raw data, or participant identifiers.

set.seed(42)

required_pkgs <- c("glmnet", "limma", "pROC")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(glmnet)
  library(limma)
  library(pROC)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!startsWith(args[[i]], "--")) stop("Unexpected argument: ", args[[i]])
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
PROTEIN_ANNOTATION_FILE <- get_arg("protein-annotation", Sys.getenv("PROTEIN_ANNOTATION_FILE", unset = NA))
OUT_DIR <- get_arg("out", Sys.getenv("OUT_DIR", unset = "results/lasso_multiclass"))
GROUP_LEVELS <- strsplit(get_arg("classes", Sys.getenv("GROUP_LEVELS", unset = "NC,HM,UM")), ",", fixed = TRUE)[[1]]
GROUP_LEVELS <- trimws(GROUP_LEVELS)

N_REPEATS <- as.integer(get_arg("repeats", Sys.getenv("N_REPEATS", unset = "10")))
K_FOLD <- as.integer(get_arg("folds", Sys.getenv("K_FOLD", unset = "5")))
ALPHA <- as.numeric(get_arg("alpha", Sys.getenv("ALPHA", unset = "1")))
PREFILTER_FDR <- as.numeric(get_arg("prefilter-fdr", Sys.getenv("PREFILTER_FDR", unset = "0.05")))
PREFILTER_FALLBACK_P <- as.numeric(get_arg("prefilter-fallback-p", Sys.getenv("PREFILTER_FALLBACK_P", unset = "0.01")))
MAX_CANDIDATES_PER_OMICS <- as.integer(get_arg("max-candidates", Sys.getenv("MAX_CANDIDATES", unset = "120")))
MIN_FINAL_FEATURES <- as.integer(get_arg("min-final-features", Sys.getenv("MIN_FINAL_FEATURES", unset = "3")))

if (is.na(PHENOTYPE_FILE) || !file.exists(PHENOTYPE_FILE)) {
  stop("Provide a phenotype CSV with --phenotype or PHENOTYPE_FILE.")
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_feature_matrix <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) < 2) stop("Feature matrix must contain feature IDs plus sample columns: ", path)
  feature_id <- make.unique(as.character(df[[1]]))
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- feature_id
  mat
}

phenotype <- read.csv(PHENOTYPE_FILE, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample_id", "class_label") %in% names(phenotype))) {
  stop("Phenotype CSV must contain columns: sample_id, class_label.")
}
phenotype$sample_id <- as.character(phenotype$sample_id)
phenotype$class_label <- factor(as.character(phenotype$class_label), levels = GROUP_LEVELS)
if (anyNA(phenotype$class_label)) {
  stop("Some class_label values are missing from --classes: ",
       paste(unique(phenotype$class_label[is.na(phenotype$class_label)]), collapse = ", "))
}
if (anyDuplicated(phenotype$sample_id)) stop("sample_id values must be unique.")

omics_mats <- list(
  Metabolomics = read_feature_matrix(METABOLOMICS_FILE),
  Lipidomics = read_feature_matrix(LIPIDOMICS_FILE),
  Proteomics = read_feature_matrix(PROTEOMICS_FILE)
)
omics_mats <- omics_mats[!vapply(omics_mats, is.null, logical(1))]
if (length(omics_mats) == 0) stop("Provide at least one omics matrix.")

for (omics_name in names(omics_mats)) {
  missing_samples <- setdiff(phenotype$sample_id, colnames(omics_mats[[omics_name]]))
  if (length(missing_samples) > 0) {
    stop(omics_name, " matrix is missing sample IDs, e.g. ", paste(head(missing_samples, 5), collapse = ", "))
  }
  omics_mats[[omics_name]] <- omics_mats[[omics_name]][, phenotype$sample_id, drop = FALSE]
}

feature_annotation <- NULL
if (!is.na(PROTEIN_ANNOTATION_FILE) && file.exists(PROTEIN_ANNOTATION_FILE)) {
  feature_annotation <- read.csv(PROTEIN_ANNOTATION_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("feature_id", "display_name") %in% names(feature_annotation))) {
    warning("Protein annotation ignored because required columns feature_id/display_name were not found.")
    feature_annotation <- NULL
  }
}

model_specs <- as.list(names(omics_mats))
names(model_specs) <- names(omics_mats)
if (length(omics_mats) > 1) {
  model_specs$Multiomics <- names(omics_mats)
}

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  y <- as.character(y)
  fold <- rep(NA_integer_, length(y))
  for (level in unique(y)) {
    idx <- sample(which(y == level))
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

prefilter_candidates <- function(mat, y) {
  y <- factor(y, levels = GROUP_LEVELS)
  if (length(unique(y)) < 2) stop("At least two classes are required in the training data.")
  design <- model.matrix(~ y)
  fit <- limma::eBayes(limma::lmFit(as.matrix(mat), design))
  res <- limma::topTable(fit, coef = 2:ncol(design), number = Inf,
                         adjust.method = "BH", sort.by = "F")
  res <- res[order(res$P.Value), , drop = FALSE]
  selected <- rownames(res)[res$adj.P.Val < PREFILTER_FDR]
  if (length(selected) < MIN_FINAL_FEATURES) {
    selected <- rownames(res)[res$P.Value < PREFILTER_FALLBACK_P]
  }
  if (length(selected) == 0) {
    selected <- rownames(res)[seq_len(min(MAX_CANDIDATES_PER_OMICS, nrow(res)))]
  }
  head(selected, MAX_CANDIDATES_PER_OMICS)
}

build_raw_matrix <- function(sample_idx, selected, model_omics) {
  parts <- list()
  for (omics_name in model_omics) {
    features <- selected[[omics_name]]
    if (length(features) == 0) next
    mat <- omics_mats[[omics_name]][features, sample_idx, drop = FALSE]
    x_part <- t(mat)
    colnames(x_part) <- paste(omics_name, features, sep = "__")
    parts[[omics_name]] <- x_part
  }
  if (length(parts) == 0) return(NULL)
  do.call(cbind, parts)
}

fit_preprocessor <- function(x_train) {
  x_train <- as.matrix(x_train)
  center <- colMeans(x_train, na.rm = TRUE)
  center[!is.finite(center)] <- 0
  imputed <- x_train
  for (j in seq_len(ncol(imputed))) {
    missing <- is.na(imputed[, j])
    if (any(missing)) imputed[missing, j] <- center[[j]]
  }
  scale <- apply(imputed, 2, sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(center = center, scale = scale)
}

apply_preprocessor <- function(x, pp) {
  x <- as.matrix(x)
  for (j in seq_len(ncol(x))) {
    missing <- is.na(x[, j])
    if (any(missing)) x[missing, j] <- pp$center[[j]]
  }
  sweep(sweep(x, 2, pp$center, "-"), 2, pp$scale, "/")
}

calc_balanced_accuracy <- function(true, pred) {
  recalls <- vapply(GROUP_LEVELS, function(cls) {
    idx <- true == cls
    if (!any(idx)) return(NA_real_)
    mean(pred[idx] == true[idx])
  }, numeric(1))
  mean(recalls, na.rm = TRUE)
}

calc_auc_ovr <- function(true, prob_mat) {
  aucs <- setNames(rep(NA_real_, length(GROUP_LEVELS)), GROUP_LEVELS)
  for (cls in GROUP_LEVELS) {
    y_binary <- as.integer(true == cls)
    if (length(unique(y_binary)) < 2) next
    roc_obj <- pROC::roc(y_binary, prob_mat[, cls], levels = c(0, 1),
                         direction = "<", quiet = TRUE)
    aucs[[cls]] <- as.numeric(pROC::auc(roc_obj))
  }
  aucs
}

nonzero_multinomial_features <- function(coef_list) {
  features <- Reduce(union, lapply(coef_list, function(coef_mat) {
    values <- as.vector(coef_mat)
    names <- rownames(coef_mat)
    names[values != 0 & names != "(Intercept)"]
  }))
  sort(features)
}

fit_cv_glmnet <- function(x, y, seed) {
  foldid <- make_stratified_folds(y, k = K_FOLD, seed = seed)
  glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "multinomial",
    alpha = ALPHA,
    foldid = foldid,
    type.measure = "deviance",
    standardize = FALSE
  )
}

make_sample_level_predictions <- function(pred_df) {
  prob_cols <- paste0("prob_", GROUP_LEVELS)
  agg <- aggregate(pred_df[, prob_cols, drop = FALSE],
                   by = list(sample_id = pred_df$sample_id, true_label = pred_df$true_label),
                   FUN = function(x) mean(x, na.rm = TRUE))
  names(agg)[match(prob_cols, names(agg))] <- paste0("mean_", prob_cols)
  counts <- aggregate(pred_label ~ sample_id + true_label, data = pred_df, FUN = length)
  names(counts)[names(counts) == "pred_label"] <- "n_outer_predictions"
  agg <- merge(agg, counts, by = c("sample_id", "true_label"), all.x = TRUE)
  prob_mat <- as.matrix(agg[, paste0("mean_prob_", GROUP_LEVELS), drop = FALSE])
  colnames(prob_mat) <- GROUP_LEVELS
  agg$pred_label <- GROUP_LEVELS[max.col(prob_mat, ties.method = "first")]
  agg[order(match(agg$true_label, GROUP_LEVELS), agg$sample_id), ]
}

feature_label_table <- function(raw_features) {
  split <- strsplit(raw_features, "__", fixed = TRUE)
  omics <- vapply(split, `[`, character(1), 1)
  feature_id <- vapply(split, function(x) paste(x[-1], collapse = "__"), character(1))
  display_name <- feature_id
  if (!is.null(feature_annotation)) {
    mapped <- feature_annotation$display_name[match(feature_id, feature_annotation$feature_id)]
    display_name <- ifelse(!is.na(mapped) & mapped != "", mapped, display_name)
  }
  data.frame(raw_feature = raw_features, omics = omics, feature_id = feature_id,
             display_name = display_name, stringsAsFactors = FALSE)
}

run_model <- function(model_name, model_omics) {
  message("Running model: ", model_name)
  model_dir <- file.path(OUT_DIR, model_name)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  n <- nrow(phenotype)
  y_all <- phenotype$class_label
  pred_rows <- list()
  perf_rows <- list()
  stability_features <- list()
  prefilter_rows <- list()

  for (repeat_id in seq_len(N_REPEATS)) {
    outer_fold <- make_stratified_folds(y_all, k = K_FOLD, seed = 1000 + repeat_id)
    repeat_true <- rep(NA_character_, n)
    repeat_pred <- rep(NA_character_, n)
    repeat_prob <- matrix(NA_real_, nrow = n, ncol = length(GROUP_LEVELS),
                          dimnames = list(NULL, GROUP_LEVELS))

    for (fold in seq_len(K_FOLD)) {
      test_idx <- which(outer_fold == fold)
      train_idx <- which(outer_fold != fold)
      y_train <- y_all[train_idx]

      selected <- setNames(vector("list", length(model_omics)), model_omics)
      for (omics_name in model_omics) {
        selected[[omics_name]] <- prefilter_candidates(omics_mats[[omics_name]][, train_idx, drop = FALSE], y_train)
      }
      prefilter_rows[[length(prefilter_rows) + 1]] <- data.frame(
        model = model_name,
        repeat_id = repeat_id,
        fold = fold,
        omics = names(selected),
        n_candidates = vapply(selected, length, integer(1)),
        stringsAsFactors = FALSE
      )

      x_train_raw <- build_raw_matrix(train_idx, selected, model_omics)
      x_test_raw <- build_raw_matrix(test_idx, selected, model_omics)
      pp <- fit_preprocessor(x_train_raw)
      x_train <- apply_preprocessor(x_train_raw, pp)
      x_test <- apply_preprocessor(x_test_raw, pp)

      cvfit <- fit_cv_glmnet(x_train, y_train, seed = 10000 + repeat_id * 100 + fold)
      prob <- predict(cvfit, newx = x_test, s = cvfit$lambda.1se, type = "response")[, , 1]
      prob <- prob[, GROUP_LEVELS, drop = FALSE]
      pred <- colnames(prob)[max.col(prob, ties.method = "first")]

      repeat_true[test_idx] <- as.character(y_all[test_idx])
      repeat_pred[test_idx] <- pred
      repeat_prob[test_idx, ] <- prob
      stability_features[[length(stability_features) + 1]] <- nonzero_multinomial_features(coef(cvfit, s = "lambda.1se"))

      pred_rows[[length(pred_rows) + 1]] <- data.frame(
        model = model_name,
        repeat_id = repeat_id,
        fold = fold,
        sample_id = phenotype$sample_id[test_idx],
        true_label = as.character(y_all[test_idx]),
        pred_label = pred,
        prob,
        lambda_min = cvfit$lambda.min,
        lambda_1se = cvfit$lambda.1se,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      names(pred_rows[[length(pred_rows)]])[match(GROUP_LEVELS, names(pred_rows[[length(pred_rows)]]))] <- paste0("prob_", GROUP_LEVELS)
    }

    valid <- !is.na(repeat_true)
    aucs <- calc_auc_ovr(repeat_true[valid], repeat_prob[valid, , drop = FALSE])
    perf_rows[[length(perf_rows) + 1]] <- data.frame(
      model = model_name,
      repeat_id = repeat_id,
      n_valid = sum(valid),
      accuracy = mean(repeat_pred[valid] == repeat_true[valid]),
      balanced_accuracy = calc_balanced_accuracy(repeat_true[valid], repeat_pred[valid]),
      macro_auc = mean(aucs, na.rm = TRUE),
      t(aucs),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  predictions <- do.call(rbind, pred_rows)
  performance <- do.call(rbind, perf_rows)
  prefilter_counts <- do.call(rbind, prefilter_rows)
  sample_predictions <- make_sample_level_predictions(predictions)

  write.csv(predictions, file.path(model_dir, paste0(model_name, "_nested_cv_predictions.csv")), row.names = FALSE)
  write.csv(performance, file.path(model_dir, paste0(model_name, "_performance_by_repeat.csv")), row.names = FALSE)
  write.csv(prefilter_counts, file.path(model_dir, paste0(model_name, "_fold_prefilter_candidate_counts.csv")), row.names = FALSE)
  write.csv(sample_predictions, file.path(model_dir, paste0(model_name, "_sample_level_predictions.csv")), row.names = FALSE)

  freq <- sort(table(unlist(stability_features)), decreasing = TRUE)
  if (length(freq) > 0) {
    stability <- cbind(
      feature_label_table(names(freq)),
      selected_times = as.integer(freq),
      total_outer_folds = length(stability_features),
      selection_frequency = as.integer(freq) / length(stability_features)
    )
  } else {
    stability <- data.frame()
  }
  write.csv(stability, file.path(model_dir, paste0(model_name, "_feature_stability.csv")), row.names = FALSE)

  selected_full <- setNames(vector("list", length(model_omics)), model_omics)
  for (omics_name in model_omics) {
    selected_full[[omics_name]] <- prefilter_candidates(omics_mats[[omics_name]], y_all)
  }
  x_full_raw <- build_raw_matrix(seq_len(n), selected_full, model_omics)
  pp_full <- fit_preprocessor(x_full_raw)
  x_full <- apply_preprocessor(x_full_raw, pp_full)
  final_cv <- fit_cv_glmnet(x_full, y_all, seed = 9000 + match(model_name, names(model_specs)))
  final_fit <- glmnet::glmnet(x_full, y_all, family = "multinomial", alpha = ALPHA,
                              lambda = final_cv$lambda.1se, standardize = FALSE)
  coef_list <- coef(final_fit, s = final_cv$lambda.1se)
  final_features <- nonzero_multinomial_features(coef_list)
  if (length(final_features) < MIN_FINAL_FEATURES) {
    final_fit <- glmnet::glmnet(x_full, y_all, family = "multinomial", alpha = ALPHA,
                                lambda = final_cv$lambda.min, standardize = FALSE)
    coef_list <- coef(final_fit, s = final_cv$lambda.min)
    final_features <- nonzero_multinomial_features(coef_list)
  }

  coef_matrix <- matrix(0, nrow = length(final_features), ncol = length(GROUP_LEVELS),
                        dimnames = list(final_features, GROUP_LEVELS))
  for (cls in GROUP_LEVELS) {
    coef_matrix[, cls] <- as.vector(coef_list[[cls]][final_features, 1])
  }
  final_panel <- cbind(feature_label_table(final_features), as.data.frame(coef_matrix, check.names = FALSE))
  write.csv(final_panel, file.path(model_dir, paste0(model_name, "_final_feature_panel.csv")), row.names = FALSE)

  list(predictions = predictions, performance = performance, sample_predictions = sample_predictions)
}

all_results <- lapply(names(model_specs), function(model_name) run_model(model_name, model_specs[[model_name]]))
names(all_results) <- names(model_specs)

all_performance <- do.call(rbind, lapply(all_results, `[[`, "performance"))
write.csv(all_performance, file.path(OUT_DIR, "all_models_performance_by_repeat.csv"), row.names = FALSE)

metric_summary <- do.call(rbind, lapply(split(all_performance, all_performance$model), function(df) {
  metrics <- setdiff(names(df), c("model", "repeat_id", "n_valid"))
  do.call(rbind, lapply(metrics, function(metric) {
    values <- df[[metric]]
    data.frame(
      model = df$model[[1]],
      metric = metric,
      mean = mean(values, na.rm = TRUE),
      sd = sd(values, na.rm = TRUE),
      q025 = unname(quantile(values, 0.025, na.rm = TRUE)),
      q975 = unname(quantile(values, 0.975, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(metric_summary, file.path(OUT_DIR, "all_models_performance_summary.csv"), row.names = FALSE)

message("Completed multiclass LASSO workflow. Outputs written to: ", normalizePath(OUT_DIR, mustWork = FALSE))
