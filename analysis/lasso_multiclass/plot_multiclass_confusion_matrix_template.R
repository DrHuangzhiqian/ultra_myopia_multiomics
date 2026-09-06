## Plot a sample-level multiclass confusion matrix from LASSO prediction output.

required_pkgs <- c("ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages(library(ggplot2))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!startsWith(args[[i]], "--")) stop("Unexpected argument: ", args[[i]])
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

args <- parse_args()
predictions_file <- args[["predictions"]]
out_dir <- if (is.null(args[["out"]])) "results/confusion_matrix" else args[["out"]]
classes <- strsplit(if (is.null(args[["classes"]])) "NC,HM,UM" else args[["classes"]], ",", fixed = TRUE)[[1]]
classes <- trimws(classes)

if (is.null(predictions_file) || !file.exists(predictions_file)) {
  stop("Provide sample-level predictions with --predictions.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pred <- read.csv(predictions_file, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("true_label", "pred_label")
if (!all(required_cols %in% names(pred))) {
  stop("Prediction file must contain: true_label, pred_label.")
}

cm <- as.data.frame(table(
  True = factor(pred$true_label, levels = classes),
  Predicted = factor(pred$pred_label, levels = classes)
))
names(cm)[3] <- "count"
cm$row_percent <- ave(cm$count, cm$True, FUN = function(x) if (sum(x) == 0) NA_real_ else x / sum(x))
cm$label <- sprintf("%d\n(%.0f%%)", cm$count, 100 * cm$row_percent)

accuracy <- mean(pred$true_label == pred$pred_label)
fig <- ggplot(cm, aes(x = Predicted, y = True, fill = row_percent)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3.6, fontface = "bold", lineheight = 0.95) +
  scale_fill_gradient(name = "Recall", low = "#EFF6FF", high = "#08519C",
                      limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")) +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_y_discrete(limits = rev(classes), drop = FALSE) +
  coord_fixed() +
  labs(title = "Confusion Matrix",
       subtitle = sprintf("Accuracy: %.1f%% | n = %d", 100 * accuracy, nrow(pred)),
       x = "Predicted Class", y = "True Class") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold", color = "grey20"),
        legend.title = element_text(face = "bold"))

write.csv(cm, file.path(out_dir, "sample_level_confusion_matrix.csv"), row.names = FALSE)
ggsave(file.path(out_dir, "sample_level_confusion_matrix.png"), fig, width = 4.2, height = 3.5, dpi = 300)
ggsave(file.path(out_dir, "sample_level_confusion_matrix.pdf"), fig, width = 4.2, height = 3.5)
