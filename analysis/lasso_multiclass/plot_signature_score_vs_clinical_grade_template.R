## Plot association between a model-derived class probability and an ordered clinical grade.

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
clinical_file <- args[["clinical"]]
target_class <- if (is.null(args[["target-class"]])) "UM" else args[["target-class"]]
grade_column <- if (is.null(args[["grade-column"]])) "clinical_grade" else args[["grade-column"]]
out_dir <- if (is.null(args[["out"]])) "results/clinical_anchor" else args[["out"]]

if (is.null(predictions_file) || !file.exists(predictions_file)) stop("Provide --predictions.")
if (is.null(clinical_file) || !file.exists(clinical_file)) stop("Provide --clinical.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pred <- read.csv(predictions_file, stringsAsFactors = FALSE, check.names = FALSE)
clinical <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
prob_col <- paste0("mean_prob_", target_class)

if (!all(c("sample_id", "true_label", prob_col) %in% names(pred))) {
  stop("Prediction file must contain sample_id, true_label, and ", prob_col, ".")
}
if (!all(c("sample_id", grade_column) %in% names(clinical))) {
  stop("Clinical file must contain sample_id and ", grade_column, ".")
}

merged <- merge(
  pred[, c("sample_id", "true_label", prob_col)],
  clinical[, c("sample_id", grade_column)],
  by = "sample_id",
  all = FALSE
)
names(merged)[names(merged) == prob_col] <- "signature_score"
names(merged)[names(merged) == grade_column] <- "clinical_grade"
merged$clinical_grade <- as.numeric(merged$clinical_grade)
merged <- merged[is.finite(merged$signature_score) & is.finite(merged$clinical_grade), ]

spearman <- suppressWarnings(cor.test(merged$signature_score, merged$clinical_grade,
                                      method = "spearman", exact = FALSE))
fit <- lm(clinical_grade ~ signature_score, data = merged)
fit_summary <- summary(fit)

summary_df <- data.frame(
  n = nrow(merged),
  target_class = target_class,
  spearman_rho = unname(spearman$estimate),
  spearman_p = spearman$p.value,
  linear_slope = unname(coef(fit)[2]),
  linear_p = coef(fit_summary)[2, "Pr(>|t|)"],
  linear_r_squared = fit_summary$r.squared,
  stringsAsFactors = FALSE
)
write.csv(merged, file.path(out_dir, "signature_score_clinical_grade_merged.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "signature_score_clinical_grade_summary.csv"), row.names = FALSE)

p_label <- if (spearman$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", spearman$p.value)
annotation <- sprintf("Spearman rho = %.2f, %s", unname(spearman$estimate), p_label)

fig <- ggplot(merged, aes(x = signature_score, y = clinical_grade, color = true_label)) +
  geom_jitter(width = 0, height = 0.045, size = 2.4, alpha = 0.82) +
  geom_smooth(inherit.aes = FALSE, aes(x = signature_score, y = clinical_grade),
              method = "lm", se = TRUE, color = "black", fill = "grey70",
              linewidth = 0.7, alpha = 0.22) +
  annotate("text", x = -Inf, y = Inf, label = annotation, hjust = -0.05, vjust = 1.3, size = 3.0) +
  labs(x = paste(target_class, "signature score"), y = "Clinical grade", color = "Class") +
  theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "right")

ggsave(file.path(out_dir, "signature_score_vs_clinical_grade.png"), fig, width = 4.2, height = 3.6, dpi = 300)
ggsave(file.path(out_dir, "signature_score_vs_clinical_grade.pdf"), fig, width = 4.2, height = 3.6)
