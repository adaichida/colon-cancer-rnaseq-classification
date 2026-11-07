# Colon cancer RNA-seq classification
# Run from the repository root: Rscript scripts/colon_cancer_classification.R

source("scripts/setup.R")

# Load required packages
library(SummarizedExperiment)
library(edgeR)
library(DESeq2)
library(caret)
library(C50)
library(randomForest)
library(ggplot2)



# IMPORTING AND SUBSETTING
# Import sample table
sampleTable <- read.table("data/colon_cancer_sample_table.txt",
                         header = TRUE, 
                         as.is = TRUE, 
                         row.names = 1, 
                         sep = "\t")


# Subset to FF samples only
sampleTable_FF <- sampleTable[sampleTable$Storage == "FF", ]


# Check the subset
head(sampleTable_FF)


# Define class vector (Disease column)
disease <- factor(sampleTable_FF$Disease)
levels(disease)
class(disease)



# Import count table
countTable <- read.table("data/colon_cancer_gene_counts.txt",
                         header = TRUE, 
                         as.is = TRUE, 
                         row.names = 1, 
                         sep = "\t")


# Check the sample names of countTable
head(colnames(countTable))

# Check the row names of sample table FF
rownames(sampleTable_FF)

# Subset countTable to only include FF samples
countTable_FF <- countTable[, rownames(sampleTable_FF)]

# Verify the dimensions
dim(countTable_FF)

# Verify the order matches
stopifnot(identical(colnames(countTable_FF), rownames(sampleTable_FF)))




# FILTERING LOW EXPRESSED GENES
# Save plots from this run in a single PDF.
pdf("results/analysis_plots.pdf", width = 9, height = 7)

# Calculate meanLog2CPM
meanLog2CPM <- rowMeans(log2(cpm(countTable_FF) + 1))


# Plot filtering histogram
hist(meanLog2CPM,
     breaks = 50,
     col = "green",
     main = "Distribution of meanlog2CPM values",
     xlab = "Mean log2CPM")

abline(v = 1, col = "red", lwd = 2, lty = 2)


# Keep genes with meanLog2CPM > 1
countTable_FF <- countTable_FF[meanLog2CPM > 1, ]

# Check dimensions after filtering
dim(countTable_FF)



# NORMALIZATION
# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = countTable_FF,
  colData = data.frame(
    Disease = disease, row.names = rownames(sampleTable_FF)),
  design = ~ Disease)


# Apply VST normalization
normCounts <- vst(dds, blind = FALSE)

# Check normalized counts
dim(assay(normCounts))



# SETUP MODEL FITTING FUNCTIONS for RANDOM FOREST AND DECISION TREE
# Random forest fit function
# New score function
rfSBF <- caret::rfSBF
rfSBF$score <- function(x, y){
  sd(x) / mean(x)
}


# New filter function
rfSBF$filter <- function(score, x, y){
  topCoV <- score > quantile(score, 0.9)
  names(score) %in% names(score)[topCoV]
}


# Decision tree fit function using C5.0
dtSBF <- rfSBF
dtSBF$fit <- function(x, y, ...){
  loadNamespace("C50")
  C50::C5.0(x = x, y = y)
}



# DECISION TREE
# Create training set
trainingSet <- data.frame(t(assay(normCounts)))

# Leave one sample out, including when another sample comes from the same individual.
# Filtering and VST above use the full dataset before resampling.
set.seed(123)

# Train and validate Decision Tree (C5.0)
dtControl <- sbfControl(
  functions = dtSBF,
  method = "loocv",
  saveDetails = TRUE,
  verbose = FALSE)


dtModel <- sbf(
  trainingSet,
  disease,
  sbfControl = dtControl)

# Show model
summary(dtModel$fit)


# Print confusion matrix from validation
confusionMatrix(dtModel$pred$predictions$pred,
                dtModel$pred$predictions$obs)


# OR

# Cross-validated predictions
dt_pred <- dtModel$pred$predictions$pred

# Confusion matrix
dt_cm <- confusionMatrix(dt_pred, dtModel$pred$predictions$obs)
dt_cm



# RANDOM FOREST 
# Set random seed
set.seed(123)

# Train and validate Random Forest
rfControl <- sbfControl(
  functions = rfSBF,
  method = "loocv",
  saveDetails = TRUE,
  verbose = FALSE)

rfModel <- sbf(
  trainingSet,
  disease,
  sbfControl = rfControl)

# Show model
rfModel$fit


# Print confusion matrix from validation
confusionMatrix(rfModel$pred$predictions$pred,
                rfModel$pred$predictions$obs)


# OR

# Cross-validated predictions
rf_pred <- rfModel$pred$predictions$pred

# Confusion matrix
rf_cm <- confusionMatrix(rf_pred, rfModel$pred$predictions$obs)
rf_cm


# CONFUSION MATRIX FOR PRESENTATION
plot_confusion_heatmap <- function(cm, title){
  
  cm_df <- as.data.frame(cm$table)
  
  ggplot(cm_df, aes(x = Prediction, y = Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), size = 7) +
    theme_minimal(base_size = 16) +
    labs(
      title = title,
      x = "Predicted class",
      y = "Actual class"
    )
}

print(plot_confusion_heatmap(dt_cm, "Decision tree"))

print(plot_confusion_heatmap(rf_cm, "Random forest"))



# MODEL ACCURACY 
model_comparison <- data.frame(
  Model = c("Decision Tree", "Random Forest"),
  Accuracy = c(dt_cm$overall["Accuracy"], rf_cm$overall["Accuracy"]),
  Kappa = c(dt_cm$overall["Kappa"], rf_cm$overall["Kappa"])
)

model_comparison


# PLOTTING 
accuracy_plot <- ggplot(model_comparison, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_col(width = 0.6) +
  ylim(0, 1) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none") +
  labs(
    title = "Cross-validation Accuracy Comparison",
    x = "Classification Algorithms",
    y = "Accuracy"
  )

print(accuracy_plot)



# VARIABLE IMPORTANCE
# Load visualization libraries
library(ggplot2)
library(dplyr)

# 1. Extract feature importance from your specific fitted model
importance_matrix <- as.data.frame(rfModel$fit$importance)

# 2. Format data for ggplot (assuming MeanDecreaseGini is the metric)
importance_df <- data.frame(
  Feature = rownames(importance_matrix),
  Importance = importance_matrix$MeanDecreaseGini
) %>%
  arrange(desc(Importance)) %>%
  slice_head(n = 10) # Keep the top 10 variables for presentation clarity

# 3. Create a polished, professional horizontal bar chart
importance_plot <- ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#2c3e50", width = 0.7) +
  coord_flip() +
  labs(
    x = "Features",
    y = "Variable Importance Score"
  ) +
  theme_minimal(base_size = 14) + # Clean theme with larger text for slides
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "#1a252f"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


print(importance_plot)
dev.off()

write.csv(model_comparison, "results/model_comparison.csv", row.names = FALSE)
write.csv(as.data.frame(dt_cm$table), "results/decision_tree_confusion_matrix.csv", row.names = FALSE)
write.csv(as.data.frame(rf_cm$table), "results/random_forest_confusion_matrix.csv", row.names = FALSE)
write.csv(importance_df, "results/random_forest_top_features.csv", row.names = FALSE)
saveRDS(list(decision_tree = dtModel, random_forest = rfModel), "results/models.rds")
capture.output(sessionInfo(), file = "results/session_info.txt")
