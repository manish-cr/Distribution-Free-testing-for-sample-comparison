setwd("/home/manish/Documents/NUS stuff/y4s1/DSA4288S-Mult/code/analysis")
library(tidyverse)
source("../test_gfs.R")
library(Boruta)
# there are 357 samples and 22278 genes. 2 classes
df <- read.csv("/home/manish/Documents/NUS stuff/y4s1/DSA4288S/data/liver_cancer.csv",
               stringsAsFactors = T)

# these will be the types
table(df['type'])
str(df[,1:3])

# train-test split
set.seed(101)
size <- floor(.8*nrow(df))
sample_id <- sample.int(n=nrow(df),size=size)

train <- df[sample_id,]
test <- df[-sample_id,]

# doing boruta algorithm to get the estimated important features
boruta.train <- boruta.train <- Boruta(type~.-samples, data = train, doTrace = 2)

# taking just the important attributes as the signal set:
final.boruta <- getSelectedAttributes(boruta.train, withTentative = F)
print(final.boruta)

# applying this on df
train_boruta <- train[,c("type", final.boruta)]
test_boruta <- test[,c("type", final.boruta)]

# Using random forest (baseline)
# Load the randomForest package
library(randomForest)

# Step 1: Fit Random Forest model on the training set
rf_model_boruta <- randomForest(type ~ ., data = train_boruta, ntree = 500, mtry = sqrt(ncol(train_boruta) - 1), importance = TRUE)

# Step 2: Predict on the test set
rf_predicted_classes <- predict(rf_model_boruta, newdata = test_boruta)

# Step 3: Calculate the test error
actual_classes <- test_boruta$type  # True labels in the test set
rf_test_error_boruta <- mean(rf_predicted_classes != actual_classes)  # Misclassification error rate

# Print the Random Forest test error
print(paste("Random Forest Test error:", rf_test_error_boruta))

# Optional: Check variable importance
importance(rf_model_boruta)

#######################################
# Comparing with features yielded from GFS
HCC_indices <- which(df$type=="HCC")
normal_indices <- which(train$type == "normal")

# PROBLEM WITH MCM
# Extract the features (excluding samples and type)
data_matrix_HCC <- as.matrix(train[HCC_indices, -(1:2)])  # Exclude the first two columns (type and samples)
data_matrix_normal <- as.matrix(train[normal_indices, -(1:2)])

# set the data list
data_list <- list(data_matrix_HCC,
                  data_matrix_normal)

# Feature set: All gene indices
features <- 1:(ncol(train)-2)  # Exclude samples and type columns

# Apply GFS
alpha_thresh <- 0.05
# selected_features_mcm <- GFS(features, alpha_thresh, data_list, test_type = "MCM")
# selected_features_mmcm <- GFS(features, alpha_thresh, data_list, test_type = "MMCM")


# --------------parallelization test--------------------------------------------

# Parallelization test
library(parallel)

# Wrapper function to call GFS
gfs_wrapper <- function(feature_set) {
  GFS(features = feature_set, alpha_thresh = alpha_thresh, data_list = data_list, test_type = "MMCM")
}

# Parallel feature selection using mclapply
selected_features_list <- mclapply(list(features), gfs_wrapper, mc.cores = detectCores() - 1)

# View the selected features
print(selected_features_list)

#########################################################
# validation error check (MMCM)
# Convert numeric indices to column names for both train and test data
corrected_features <- selected_features_list[[1]] + 2 # (to adjust for the type and samples)
selected_feature_names <- c("type", colnames(train)[corrected_features])

# Subset train and test datasets using the selected feature names
train_mmcm <- train[, selected_feature_names, drop = FALSE]
test_mmcm <- test[, selected_feature_names, drop = FALSE]

# Train the Random Forest model on the selected features
rf_model_mmcm <- randomForest(type ~ ., data = train_mmcm, ntree = 500, mtry = sqrt(ncol(train_mmcm) - 1), importance = TRUE)

rf_predicted_classes_mmcm <- predict(rf_model_mmcm, newdata = test_mmcm)

actual_classes_mmcm <- test_mmcm$type  # True labels in the test set
rf_test_error_mmcm <- mean(rf_predicted_classes_mmcm != actual_classes_mmcm)  # Misclassification error rate

# Print the Random Forest test error
print(paste("Random Forest Test error:", rf_test_error_mmcm))
