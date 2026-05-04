##---------------------
# Improved Linear Discriminant Analysis for Stomata Ploidy Prediction
# Author: Improved version
# Date: July 15, 2025
##---------------------

# Enable the r-universe repo
options(repos = c(
  fawda123 = 'https://fawda123.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

# Load required libraries
library(readr)
library(klaR)
library(psych)
library(ggord)
library(devtools)
library(MASS)
library(ggplot2)
library(caret)
library(dplyr)
library(gridExtra)

# Set working directory (adjust as needed)
setwd("C:/Users/Jaume/Nextcloud/Documents/Pol/Tmesepteris")

##---------------------
# DATA LOADING AND INITIAL EXPLORATION
##---------------------

# Load training data (with known ploidy levels)
train_df <- read_csv("C:/Users/Jaume/Nextcloud/Documents/Pol/Tmesepteris/stomata_measured_ploidy.csv", )

# Load prediction data (unknown ploidy levels to be predicted)
predict_df <- read_csv("C:/Users/Jaume/Nextcloud/Documents/Pol/Tmesepteris/data_model_predictive.csv")
# Remove columns with all NA values from training data
na_cols_train <- sapply(train_df, function(x) all(is.na(x)))
if(any(na_cols_train)) {
  cat("Removing columns with all NA values from training data:", names(train_df)[na_cols_train], "\n")
  train_df <- train_df[, !na_cols_train]
}

# Remove columns with all NA values from prediction data
na_cols_pred <- sapply(predict_df, function(x) all(is.na(x)))
if(any(na_cols_pred)) {
  cat("Removing columns with all NA values from prediction data:", names(predict_df)[na_cols_pred], "\n")
  predict_df <- predict_df[, !na_cols_pred]
}
# Display basic information about the training dataset
cat("=== TRAINING DATASET ===\n")
cat("Dataset dimensions:", dim(train_df), "\n")
cat("Column names:", names(train_df), "\n")
cat("Summary statistics:\n")
summary(train_df)

# Check for missing values in training data
cat("Missing values per column in training data:\n")
print(colSums(is.na(train_df)))

# Display basic information about the prediction dataset
cat("\n=== PREDICTION DATASET ===\n")
cat("Dataset dimensions:", dim(predict_df), "\n")
cat("Column names:", names(predict_df), "\n")
cat("Summary statistics:\n")
summary(predict_df)

# Check for missing values in prediction data
cat("Missing values per column in prediction data:\n")
print(colSums(is.na(predict_df)))

##---------------------
# OUTLIER DETECTION AND REMOVAL
##---------------------

# Function to detect extreme outliers using IQR method (only very extreme ones)
detect_extreme_outliers <- function(x, variable_name, dataset_name) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  
  # Use 3*IQR instead of 1.5*IQR to catch only really extreme outliers
  lower_bound <- Q1 - 3 * IQR
  upper_bound <- Q3 + 3 * IQR
  
  extreme_outliers <- x[x < lower_bound | x > upper_bound]
  cat("Number of extreme outliers in", variable_name, "in", dataset_name, ":", length(extreme_outliers), "\n")
  cat("Range of", variable_name, ":", round(min(x, na.rm = TRUE), 2), "to", round(max(x, na.rm = TRUE), 2), "\n")
  cat("Extreme outlier bounds: <", round(lower_bound, 2), "or >", round(upper_bound, 2), "\n")
  
  return(extreme_outliers)
}

# Process training data
cat("=== PROCESSING TRAINING DATA ===\n")
# Detect extreme outliers in training data
train_width_outliers <- detect_extreme_outliers(train_df$Width, "Width", "training data")
train_length_outliers <- detect_extreme_outliers(train_df$Length, "Length", "training data")

# Create a copy of the training data for cleaning
train_clean <- train_df

# Remove extreme outliers from training data (only if they exist)
if(length(train_width_outliers) > 0) {
  train_clean <- train_clean[-which(train_clean$Width %in% train_width_outliers), ]
}
if(length(train_length_outliers) > 0) {
  train_clean <- train_clean[-which(train_clean$Length %in% train_length_outliers), ]
}

# Convert Ploidy to factor in training data
train_clean$Ploidy <- as.factor(train_clean$Ploidy)

cat("Training data dimensions after outlier removal:", dim(train_clean), "\n")
cat("Ploidy levels in training data:", levels(train_clean$Ploidy), "\n")

# Process prediction data
cat("\n=== PROCESSING PREDICTION DATA ===\n")
# Detect extreme outliers in prediction data
pred_width_outliers <- detect_extreme_outliers(predict_df$Width, "Width", "prediction data")
pred_length_outliers <- detect_extreme_outliers(predict_df$Length, "Length", "prediction data")

# Create a copy of the prediction data for cleaning
predict_clean <- predict_df

# Remove extreme outliers from prediction data (only if they exist)
if(length(pred_width_outliers) > 0) {
  predict_clean <- predict_clean[-which(predict_clean$Width %in% pred_width_outliers), ]
}
if(length(pred_length_outliers) > 0) {
  predict_clean <- predict_clean[-which(predict_clean$Length %in% pred_length_outliers), ]
}

# If there's a Ploidy column in prediction data, convert to factor (it might be all NAs)
if("Ploidy" %in% names(predict_clean)) {
  predict_clean$Ploidy <- as.factor(predict_clean$Ploidy)
}

cat("Prediction data dimensions after outlier removal:", dim(predict_clean), "\n")

##---------------------
# DATA PREPARATION FOR MODELING
##---------------------

# Use complete cases from training data for model training
complete_training_data <- train_clean[complete.cases(train_clean), ]

# Check if we have enough data for each ploidy level
cat("=== TRAINING DATA SUMMARY ===\n")
cat("Complete cases for training:", nrow(complete_training_data), "\n")
cat("Distribution of ploidy levels in training data:\n")
print(table(complete_training_data$Ploidy))

# Prepare prediction data (should have Length and Width available)
prediction_data <- predict_clean[!is.na(predict_clean$Length) & !is.na(predict_clean$Width), ]

cat("\n=== PREDICTION DATA SUMMARY ===\n")
cat("Cases available for prediction:", nrow(prediction_data), "\n")

# Check minimum sample size per group for LDA
min_samples <- min(table(complete_training_data$Ploidy))
cat("Minimum samples per ploidy level:", min_samples, "\n")
if(min_samples < 3) {
  cat("Warning: Some ploidy levels have very few samples. Consider combining groups or collecting more data.\n")
}

##---------------------
# EXPLORATORY DATA ANALYSIS
##---------------------

# Create comprehensive visualization of training data
p1 <- ggplot(data = complete_training_data, aes(x = Length, y = Width, color = Ploidy)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_classic() +
  labs(title = "Training Data: Stomata Measurements by Ploidy Level",
       x = "Length (μm)",
       y = "Width (μm)",
       color = "Ploidy") +
  theme(plot.title = element_text(hjust = 0.5))

# Box plots for each variable by ploidy in training data
p2 <- ggplot(complete_training_data, aes(x = Ploidy, y = Length, fill = Ploidy)) +
  geom_boxplot(alpha = 0.7) +
  theme_classic() +
  labs(title = "Training Data: Length Distribution by Ploidy",
       x = "Ploidy",
       y = "Length (μm)") +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(complete_training_data, aes(x = Ploidy, y = Width, fill = Ploidy)) +
  geom_boxplot(alpha = 0.7) +
  theme_classic() +
  labs(title = "Training Data: Width Distribution by Ploidy",
       x = "Ploidy",
       y = "Width (μm)") +
  theme(plot.title = element_text(hjust = 0.5))

# Scatter plot of prediction data (without ploidy colors)
p4 <- ggplot(data = prediction_data, aes(x = Length, y = Width)) +
  geom_point(size = 2, alpha = 0.7, color = "gray50") +
  theme_classic() +
  labs(title = "Prediction Data: Stomata Measurements",
       x = "Length (μm)",
       y = "Width (μm)") +
  theme(plot.title = element_text(hjust = 0.5))

# Display plots
print(p1)
print(p2)
print(p3)
print(p4)

# Print descriptive statistics by group for training data
cat("\nDescriptive statistics by Ploidy level (Training Data):\n")
print(describeBy(complete_training_data[, c("Length", "Width")], complete_training_data$Ploidy))

# Print descriptive statistics for prediction data
cat("\nDescriptive statistics for Prediction Data:\n")
print(describe(prediction_data[, c("Length", "Width")]))

##---------------------
# LINEAR DISCRIMINANT ANALYSIS
##---------------------

# Fit LDA model using training data (simple model with Length and Width)
lda_model <- lda(Ploidy ~ Length + Width, data = complete_training_data)

# Display LDA results
cat("\nLDA Model Summary:\n")
print(lda_model)

# Calculate proportions of trace
prop_trace <- lda_model$svd^2 / sum(lda_model$svd^2)
cat("\nProportion of trace:\n")
print(prop_trace)

##---------------------
# MODEL VALIDATION AND PERFORMANCE
##---------------------

# Cross-validation using caret
set.seed(123)  # For reproducibility
train_control <- trainControl(method = "cv", number = 10)

# Cross-validation on training data
lda_cv <- train(Ploidy ~ Length + Width, 
                method = "lda", 
                data = complete_training_data,
                trControl = train_control)

cat("\nCross-validation results:\n")
print(lda_cv)

# Confusion matrix for training data
predictions_train <- predict(lda_model, complete_training_data)$class
confusion_matrix <- confusionMatrix(predictions_train, complete_training_data$Ploidy)
cat("\nConfusion Matrix for Training Data:\n")
print(confusion_matrix)

##---------------------
# VISUALIZATION OF LDA RESULTS
##---------------------

# Get LDA values for plotting training data
lda_values <- predict(lda_model, complete_training_data)

# Create data frame for plotting
lda_df <- data.frame(
  Ploidy = complete_training_data$Ploidy,
  LD1 = lda_values$x[, 1],
  LD2 = if(ncol(lda_values$x) > 1) lda_values$x[, 2] else rep(0, nrow(lda_values$x))
)

# LDA plot
p5 <- ggplot(lda_df, aes(x = LD1, y = LD2, color = Ploidy)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_classic() +
  labs(title = "LDA Results - Training Data",
       x = paste("LD1 (", round(prop_trace[1] * 100, 1), "% of variance)", sep = ""),
       y = if(length(prop_trace) > 1) paste("LD2 (", round(prop_trace[2] * 100, 1), "% of variance)", sep = "") else "LD2",
       color = "Ploidy") +
  theme(plot.title = element_text(hjust = 0.5))

print(p5)

# Histogram of LD1 scores
if(length(unique(complete_training_data$Ploidy)) > 1) {
  ldahist(lda_values$x[, 1], g = complete_training_data$Ploidy, main = "LDA Histogram - Training Data")
}

# Partition plot
partimat(Ploidy ~ Length + Width, data = complete_training_data, method = "lda",
         main = "LDA Partition Plot - Training Data")

##---------------------
# PREDICTION ON NEW DATA
##---------------------

# Make predictions on the prediction dataset
predictions_new <- predict(lda_model, prediction_data)

# Add predictions to the prediction dataset
prediction_data$predicted_ploidy <- predictions_new$class
prediction_data$prediction_probability <- apply(predictions_new$posterior, 1, max)

# Get the LDA scores for the prediction data
prediction_lda_values <- predict(lda_model, prediction_data)

# Create visualization of predictions
prediction_lda_df <- data.frame(
  LD1 = prediction_lda_values$x[, 1],
  LD2 = if(ncol(prediction_lda_values$x) > 1) prediction_lda_values$x[, 2] else rep(0, nrow(prediction_lda_values$x)),
  PredictedPloidy = prediction_data$predicted_ploidy,
  Length = prediction_data$Length,
  Width = prediction_data$Width,
  Probability = prediction_data$prediction_probability
)

# Plot predictions in LDA space
p6 <- ggplot(prediction_lda_df, aes(x = LD1, y = LD2, color = PredictedPloidy)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_classic() +
  labs(title = "Predicted Ploidy Levels - LDA Space",
       x = paste("LD1 (", round(prop_trace[1] * 100, 1), "% of variance)", sep = ""),
       y = if(length(prop_trace) > 1) paste("LD2 (", round(prop_trace[2] * 100, 1), "% of variance)", sep = "") else "LD2",
       color = "Predicted Ploidy") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot predictions in original space
p7 <- ggplot(prediction_data, aes(x = Length, y = Width, color = predicted_ploidy)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_classic() +
  labs(title = "Predicted Ploidy Levels - Original Space",
       x = "Length (μm)",
       y = "Width (μm)",
       color = "Predicted Ploidy") +
  theme(plot.title = element_text(hjust = 0.5))

print(p6)
print(p7)

# Display prediction results
cat("\nPrediction Summary:\n")
prediction_table <- table(prediction_data$predicted_ploidy)
print(prediction_table)

# Show prediction confidence
cat("\nPrediction Confidence Statistics:\n")
cat("Mean prediction probability:", round(mean(prediction_data$prediction_probability), 3), "\n")
cat("Median prediction probability:", round(median(prediction_data$prediction_probability), 3), "\n")
cat("Min prediction probability:", round(min(prediction_data$prediction_probability), 3), "\n")
cat("Max prediction probability:", round(max(prediction_data$prediction_probability), 3), "\n")

# Show low-confidence predictions
low_confidence_threshold <- 0.7
low_confidence_predictions <- prediction_data[prediction_data$prediction_probability < low_confidence_threshold, ]
cat("\nLow confidence predictions (probability <", low_confidence_threshold, "):", nrow(low_confidence_predictions), "\n")

# Save predictions
write.csv(prediction_data, "predicted_stomata_ploidy.csv", row.names = FALSE)
cat("\nAll predictions saved to 'predicted_stomata_ploidy.csv'\n")

# Save low confidence predictions separately
if(nrow(low_confidence_predictions) > 0) {
  write.csv(low_confidence_predictions, "low_confidence_predictions.csv", row.names = FALSE)
  cat("Low confidence predictions saved to 'low_confidence_predictions.csv'\n")
}

##---------------------
# MODEL DIAGNOSTICS
##---------------------

# Test assumptions on training data
cat("\nTesting multivariate normality and homogeneity of variance (Training Data):\n")

# Shapiro-Wilk test for normality (for each group)
for(level in levels(complete_training_data$Ploidy)) {
  subset_data <- complete_training_data[complete_training_data$Ploidy == level, ]
  cat("\nPloidy level:", level, "\n")
  cat("Sample size:", nrow(subset_data), "\n")
  cat("Length normality test p-value:", shapiro.test(subset_data$Length)$p.value, "\n")
  cat("Width normality test p-value:", shapiro.test(subset_data$Width)$p.value, "\n")
}

# Bartlett's test for homogeneity of variance
cat("\nBartlett's test for homogeneity of variance:\n")
cat("Length:", bartlett.test(Length ~ Ploidy, data = complete_training_data)$p.value, "\n")
cat("Width:", bartlett.test(Width ~ Ploidy, data = complete_training_data)$p.value, "\n")

##---------------------
# SUMMARY AND RECOMMENDATIONS
##---------------------

cat("\n==== ANALYSIS SUMMARY ====\n")
cat("Training dataset samples:", nrow(train_df), "\n")
cat("Training samples after outlier removal:", nrow(train_clean), "\n")
cat("Complete training cases used:", nrow(complete_training_data), "\n")
cat("Prediction dataset samples:", nrow(predict_df), "\n")
cat("Prediction samples after outlier removal:", nrow(predict_clean), "\n")
cat("Cases predicted:", nrow(prediction_data), "\n")
cat("Model accuracy on training data:", round(confusion_matrix$overall['Accuracy'], 3), "\n")
cat("Model Kappa on training data:", round(confusion_matrix$overall['Kappa'], 3), "\n")
cat("Cross-validation accuracy:", round(lda_cv$results$Accuracy, 3), "\n")

##---------------------
# INDIVIDUAL-LEVEL PLOIDY ASSIGNMENT
##---------------------

cat("\n==== INDIVIDUAL-LEVEL ANALYSIS ====\n")

# Check if code variable exists
if("code" %in% names(prediction_data)) {
  
  # Calculate individual-level statistics
  individual_stats <- prediction_data %>%
    group_by(code) %>%
    summarise(
      total_stomata = n(),
      predicted_4x = sum(predicted_ploidy == "4x"),
      predicted_8x = sum(predicted_ploidy == "8x"),
      mean_probability = mean(prediction_probability),
      median_probability = median(prediction_probability),
      min_probability = min(prediction_probability),
      .groups = 'drop'
    )
  
  # Calculate percentages for each ploidy level
  individual_stats$perc_4x <- round((individual_stats$predicted_4x / individual_stats$total_stomata) * 100, 1)
  individual_stats$perc_8x <- round((individual_stats$predicted_8x / individual_stats$total_stomata) * 100, 1)
  
  # Apply 80% threshold rule for individual assignment
  individual_stats$assigned_ploidy <- NA
  individual_stats$assignment_confidence <- NA
  
  threshold <- 80
  
  for(i in 1:nrow(individual_stats)) {
    percentages <- c(individual_stats$perc_4x[i], individual_stats$perc_8x[i])
    ploidy_levels <- c("4x", "8x")
    
    max_perc <- max(percentages)
    max_ploidy <- ploidy_levels[which.max(percentages)]
    
    if(max_perc >= threshold) {
      individual_stats$assigned_ploidy[i] <- max_ploidy
      individual_stats$assignment_confidence[i] <- "High"
    } else {
      individual_stats$assigned_ploidy[i] <- "Ambiguous"
      individual_stats$assignment_confidence[i] <- "Low"
    }
  }
  
  # Display results
  cat("Individual-level ploidy assignments (≥80% threshold):\n")
  print(individual_stats[, c("code", "total_stomata", "perc_4x", "perc_8x", 
                             "assigned_ploidy", "assignment_confidence", "mean_probability")])
  
  # Summary statistics
  cat("\n==== INDIVIDUAL ASSIGNMENT SUMMARY ====\n")
  assignment_summary <- table(individual_stats$assigned_ploidy)
  print(assignment_summary)
  
  high_confidence_count <- sum(individual_stats$assignment_confidence == "High")
  low_confidence_count <- sum(individual_stats$assignment_confidence == "Low")
  
  cat("\nHigh confidence assignments:", high_confidence_count, "/", nrow(individual_stats), 
      "(", round(high_confidence_count/nrow(individual_stats)*100, 1), "%)\n")
  cat("Low confidence assignments:", low_confidence_count, "/", nrow(individual_stats), 
      "(", round(low_confidence_count/nrow(individual_stats)*100, 1), "%)\n")
  
  # Create visualization of individual assignments
  library(ggplot2)
  
  # Prepare data for plotting
  individual_plot_data <- individual_stats %>%
    select(code, perc_4x, perc_8x, assigned_ploidy, assignment_confidence) %>%
    pivot_longer(cols = c(perc_4x, perc_8x), 
                 names_to = "ploidy_level", values_to = "percentage") %>%
    mutate(ploidy_level = gsub("perc_", "", ploidy_level))
  
  # Plot individual assignments
  p8 <- ggplot(individual_plot_data, aes(x = code, y = percentage, fill = ploidy_level)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red", size = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Individual Ploidy Predictions",
         subtitle = paste("Red line shows", threshold, "% threshold for assignment"),
         x = "Individual Code",
         y = "Percentage of Stomata",
         fill = "Predicted Ploidy") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  print(p8)
  
  # Plot confidence vs number of stomata
  p9 <- ggplot(individual_stats, aes(x = total_stomata, y = mean_probability, 
                                     color = assignment_confidence, shape = assigned_ploidy)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_classic() +
    labs(title = "Individual Assignment Confidence vs Sample Size",
         x = "Number of Stomata Measured",
         y = "Mean Prediction Probability",
         color = "Assignment Confidence",
         shape = "Assigned Ploidy") +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p9)
  
  # Save individual-level results
  write.csv(individual_stats, "individual_ploidy_assignments.csv", row.names = FALSE)
  cat("\nIndividual-level assignments saved to 'individual_ploidy_assignments.csv'\n")
  
  # Save detailed individual data
  individual_detailed <- prediction_data %>%
    left_join(individual_stats[, c("code", "assigned_ploidy", "assignment_confidence")], 
              by = "code")
  
  write.csv(individual_detailed, "detailed_predictions_with_assignments.csv", row.names = FALSE)
  cat("Detailed predictions with individual assignments saved to 'detailed_predictions_with_assignments.csv'\n")
  
  # Recommendations based on results
  cat("\n==== INDIVIDUAL-LEVEL RECOMMENDATIONS ====\n")
  
  # Check sample sizes
  min_stomata <- min(individual_stats$total_stomata)
  max_stomata <- max(individual_stats$total_stomata)
  mean_stomata <- round(mean(individual_stats$total_stomata), 1)
  
  cat("Sample size per individual ranges from", min_stomata, "to", max_stomata, "stomata (mean:", mean_stomata, ")\n")
  
  if(min_stomata < 10) {
    cat("⚠️  Warning: Some individuals have very few stomata measured (< 10)\n")
    cat("   Consider measuring more stomata for better reliability\n")
  }
  
  # Check ambiguous assignments
  ambiguous_count <- sum(individual_stats$assigned_ploidy == "Ambiguous")
  if(ambiguous_count > 0) {
    cat("⚠️  ", ambiguous_count, "individuals show ambiguous ploidy patterns\n")
    cat("   These may represent:\n")
    cat("   - Mixoploid individuals (multiple ploidy levels)\n")
    cat("   - Individuals near decision boundaries\n")
    cat("   - Measurement variability\n")
  }
  
  # Check low confidence assignments
  low_conf_assignments <- individual_stats[individual_stats$assignment_confidence == "Low", ]
  if(nrow(low_conf_assignments) > 0) {
    cat("⚠️  Consider manual review of ambiguous individuals:\n")
    for(i in 1:nrow(low_conf_assignments)) {
      cat("   - Individual", low_conf_assignments$code[i], ": 4x =", 
          low_conf_assignments$perc_4x[i], "%, 8x =", 
          low_conf_assignments$perc_8x[i], "%\n")
    }
  }
  
} else {
  cat("⚠️  Warning: 'code' variable not found in prediction data\n")
  cat("   Individual-level analysis requires a 'code' column to identify individuals\n")
  cat("   Please ensure your data has a column named 'code' with individual identifiers\n")
}

#################Basic PLOT#################
library(ggplot2)
library(patchwork)
library(readxl)

# Import the Excel file
sopredata <- read_excel("C:/Users/Jaume/Nextcloud/Documents/Pol/Tmesepteris/Spore_data.xlsx")
sopredata$PL_true <- as.factor(sopredata$PL_true)


library(ggplot2)
library(patchwork)


# Base scatter plot: grey points first, then colored ones
scatter <- ggplot() +
  geom_point(data = subset(sopredata, PL_true == "na"), 
             aes(x = Length, y = Width), color = "#B3B3B3", size = 1) +
  geom_point(data = subset(sopredata, PL_true != "na"), 
             aes(x = Length, y = Width, color = PL_true), size = 1) +
  scale_color_manual(values = c("#FED976", "#B10026", "black", "#B3B3B3")) +
  theme_minimal() +
  theme(legend.position = "none")

# X-axis density plot
x_density <- ggplot(sopredata, aes(x = Length, fill = PL_true, color = PL_true)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("#FED976", "#B10026", "black", "#B3B3B3")) +
  scale_color_manual(values = c("#FED976", "#B10026", "black", "#B3B3B3")) +
  theme_void() +
  theme(legend.position = "none")

# Y-axis density plot (flipped and compressed)
y_density <- ggplot(sopredata, aes(x = Width, fill = PL_true, color = PL_true)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("#FED976", "#B10026", "black", "#B3B3B3")) +
  scale_color_manual(values = c("#FED976", "#B10026", "black", "#B3B3B3")) +
  coord_flip() +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# Combine with layout adjustments
(x_density + plot_spacer()) /
  (scatter + y_density) +
  plot_layout(heights = c(0.7, 3), widths = c(5, 0.4))
