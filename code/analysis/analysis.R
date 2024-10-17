library(tidyverse)
library(Boruta)
# there are 357 samples and 22278 genes. 2 classes
df <- read.csv("/home/manish/Documents/NUS stuff/y4s1/DSA4288S/data/liver_cancer.csv",
               stringsAsFactors = T)
# Convert specific string columns to factors
# df <- df %>% mutate_if(is.character, as.factor)

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
train <- train[,c("type", final.boruta)]
test <- test[,c("type", final.boruta)]

# Using random forest (baseline)
# Load the randomForest package
library(randomForest)

# Step 1: Fit Random Forest model on the training set
rf_model <- randomForest(type ~ ., data = train, ntree = 500, mtry = sqrt(ncol(train) - 1), importance = TRUE)

# Step 2: Predict on the test set
rf_predicted_classes <- predict(rf_model, newdata = test)

# Step 3: Calculate the test error
actual_classes <- test$type  # True labels in the test set
rf_test_error <- mean(rf_predicted_classes != actual_classes)  # Misclassification error rate

# Print the Random Forest test error
print(paste("Random Forest Test error:", rf_test_error))

# Optional: Check variable importance
importance(rf_model)

#######################################
# Comparing with features yielded from GFS
