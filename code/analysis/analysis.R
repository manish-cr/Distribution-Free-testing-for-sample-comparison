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

# Extract the features (excluding samples and type)
data_matrix_HCC <- as.matrix(train[HCC_indices, -(1:2)])  # Exclude the first two columns (type and samples)
data_matrix_normal <- as.matrix(train[normal_indices, -(1:2)])

# set the data list
data_list <- list(data_matrix_HCC,
                  data_matrix_normal)

# Feature set: All gene indices
features <- 1:(ncol(train) - 2)  # Exclude samples and type columns

# Apply GFS with "MCM"
alpha_thresh <- 0.05
selected_features_mcm <- GFS(features, alpha_thresh, data_list, test_type = "MCM")
selected_features_mmcm <- GFS(features, alpha_thresh, data_list, test_type = "MMCM")
selected_features_cf <- GFS(features, alpha_thresh, data_list, test_type = "CF")

# --------------parallelization test (FR)--------------------------------------------
# library(foreach)
# library(doParallel)
# 
# # Register the parallel backend
# num_cores <- detectCores() - 3  # Adjust the number of cores as needed
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# 
# # Define a function to process chunks of features
# process_chunk <- function(chunk_features, data_list, alpha_thresh, test_type) {
#   # Run GFS on the chunk of features
#   GFS(chunk_features, alpha_thresh, data_list, test_type)
# }
# 
# # Split features into chunks for parallel processing (e.g., 1000 features per chunk)
# chunk_size <- 1000  # Adjust the chunk size as needed
# feature_chunks <- split(features, ceiling(seq_along(features) / chunk_size))
# 
# # Use foreach to process feature chunks in parallel
# selected_features_list <- foreach(chunk = feature_chunks, .combine = c, .packages = c("multicross", "FRmatch")) %dopar% {
#   process_chunk(chunk, data_list, alpha_thresh, test_type = "FR")
# }
# 
# # Stop the cluster after processing
# stopCluster(cl)
# 
# # Print the final selected features from all chunks
# print(selected_features_list)
