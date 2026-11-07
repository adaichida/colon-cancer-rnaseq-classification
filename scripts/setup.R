# Check the packages needed by the analysis and create the output folder.
# Run this project from the folder containing README.md and setup.R.
required_packages <- c("SummarizedExperiment", "edgeR", "DESeq2", "caret",
                       "C50", "randomForest", "ggplot2", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install these packages before running the analysis: ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}
if (!file.exists("data/colon_cancer_gene_counts.txt") ||
    !file.exists("data/colon_cancer_sample_table.txt")) {
  stop("Run from the repository root, which contains the data folder.", call. = FALSE)
}
dir.create("results", showWarnings = FALSE)
