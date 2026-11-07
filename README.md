# Colon cancer RNA-seq classification

I compared decision tree and random forest models to classify carcinoma and normal colon tissue using gene expression data in R. I completed this analysis for my MSc coursework in Bioinformatics.

## Data

The supplied count matrix contains 19,351 genes across 30 samples. This analysis uses the 18 fresh-frozen samples: 12 carcinoma and 6 normal samples from 12 individuals. Some individuals contributed both tissue types.

The coursework materials link the data to Hedegaard et al. (2014), *Next-Generation Sequencing of RNA and DNA Isolated from Paired Fresh-Frozen and Formalin-Fixed Paraffin-Embedded Samples of Human Cancer and Normal Tissue*. DOI: https://doi.org/10.1371/journal.pone.0098187.

The paper lists ArrayExpress accession E-MTAB-2523 for expression matrices. This repository includes the count and sample tables supplied with the coursework. It does not recreate their upstream processing from sequencing reads.

## Analysis

1. Select fresh-frozen samples and match the sample table to the count matrix.
2. Retain genes with mean log2(CPM + 1) greater than 1.
3. Apply the DESeq2 variance-stabilising transformation.
4. Select genes above the 90th percentile of the coefficient of variation within caret's selection-by-filtering procedure.
5. Fit a C5.0 decision tree and a random forest with sample-level leave-one-out cross-validation.
6. Compare accuracy and Cohen's kappa, plot confusion matrices, and inspect random forest feature importance.

## Interpretation

The original coursework figures show 17 correct classifications out of 18 samples for each model, approximately 94.4% accuracy. The figures in `figures/` come from that original analysis. I have not rerun the adapted script to confirm those results.

The validation leaves out one sample at a time. It does not leave out an entire individual, so paired tissues can occur in different training and validation sets. Low-expression filtering and the DESeq2 transformation also use the full dataset before validation. These choices limit how strongly the accuracy estimates support prediction for new individuals. Feature importance here does not establish a validated biomarker.

## Files

- `scripts/colon_cancer_classification.R`: analysis, model comparison, and plots.
- `scripts/setup.R`: package checks and output folder creation.
- `data/`: gene counts and sample information.
- `figures/`: figures from the original coursework analysis.
- `results/`: output folder that the script creates when you run the analysis.

## Run the analysis

Keep both R scripts inside the `scripts/` folder. Run the analysis from the main repository folder, which contains `README.md`, `data/`, and `scripts/`.

Install the packages once in R:

```r
install.packages(c("BiocManager", "caret", "C50", "randomForest", "ggplot2", "dplyr"))
BiocManager::install(c("SummarizedExperiment", "edgeR", "DESeq2"))
```

Open R or RStudio in the repository folder, then run:

```r
source("scripts/colon_cancer_classification.R")
```

The analysis script loads `scripts/setup.R` to check packages and create the output folder.

Or run from a terminal in the repository folder:

```sh
Rscript scripts/colon_cancer_classification.R
```

The script writes plots, model comparisons, confusion matrices, feature importance, fitted models, and session information to `results/`.

## Code changes for this repository

The adapted script uses relative file paths, checks sample alignment, supplies plot titles, saves outputs, and draws confusion matrices from model predictions instead of manually entered counts. It retains the original sample-level validation and preprocessing order.

## Author

Jane Chidimma Ezeanowai

Msc Systems Biology with specialisation in Bioinformatics

University of Skövde
