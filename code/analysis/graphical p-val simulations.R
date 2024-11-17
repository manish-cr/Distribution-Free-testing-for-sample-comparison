setwd("/home/manish/Documents/NUS stuff/y4s1/DSA4288S-Mult/code/analysis")
source("../test_fn.R")
library(tidyverse)
library(multicross)
library(FRmatch)

# setting up baselines
# We will be using the following distributions and observing the p-value through scale, location and dimension shifts
# - normal
# - t
# - lognormal
# - weibull
# - chi squared

sep <- 5
ss <- 50

# <<<Normal dist.>>>
# location shift
# Define a function to perform the specified test for a given dimension
perform_test_norm_loc <- function(d, test_type) {
  # Generate X1 and X2 with the current dimension
  X1 <- MASS::mvrnorm(ss, rep(sep, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
  X2 <- MASS::mvrnorm(ss, rep(0, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
  
  # Perform the specified test and return the p-value
  p_val <- test_fn(test_type, list(X1, X2), 0.05)
  return(p_val)
}

# Set the range of dimensions and test types
dimensions <- 2:200
test_types <- c("MCM", "MMCM", "FR")

# Initialize a list to store p-values for each test type
p_val_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  p_val_results[[test_type]] <- mclapply(dimensions, perform_test_norm_loc, test_type = test_type, mc.cores = parallel::detectCores() - 3)
}

# Convert p-value results to a data frame for plotting
p_val_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  P_Value = unlist(p_val_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot p-values for each test type
# ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
#   geom_line() +  
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
#   labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests") +
#   theme_minimal() +
#   theme(legend.title = element_blank())

ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests using normal distribution under location shift") +
  theme_minimal() +
  theme(legend.title = element_blank())

# scale shift
perform_test_norm_scale <- function(d, test_type) {
  # Generate X1 and X2 with the current dimension
  X1 <- MASS::mvrnorm(ss, rep(0, d), diag(1+sep, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
  X2 <- MASS::mvrnorm(ss, rep(0, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
  
  # Perform the specified test and return the p-value
  p_val <- test_fn(test_type, list(X1, X2), 0.05)
  return(p_val)
}

# Set the range of dimensions and test types
dimensions <- 2:200
test_types <- c("MCM", "MMCM", "FR")

# Initialize a list to store p-values for each test type
p_val_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  p_val_results[[test_type]] <- mclapply(dimensions, perform_test_norm_scale, test_type = test_type, mc.cores = parallel::detectCores() - 3)
}

# Convert p-value results to a data frame for plotting
p_val_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  P_Value = unlist(p_val_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot p-values for each test type
# ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
#   geom_line() +  
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
#   labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests") +
#   theme_minimal() +
#   theme(legend.title = element_blank())

ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests using normal distribution under scale shift") +
  theme_minimal() +
  theme(legend.title = element_blank())
#############################################
#############################################
# <<<lognormal>>>
# location shift
# Define a function to perform the specified test for a given dimension
perform_test_lognorm_loc <- function(d, test_type) {
  # Generate X1 and X2 with the current dimension
  X1 <- MASS::mvrnorm(ss, rep(sep, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
  X2 <- MASS::mvrnorm(ss, rep(0, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
  X1 <- exp(X1)
  X2 <- exp(X2)
  # Perform the specified test and return the p-value
  p_val <- test_fn(test_type, list(X1, X2), 0.05)
  return(p_val)
}

# Set the range of dimensions and test types
dimensions <- 2:200
test_types <- c("MCM", "MMCM", "FR")

# Initialize a list to store p-values for each test type
p_val_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  p_val_results[[test_type]] <- mclapply(dimensions, perform_test_lognorm_loc, test_type = test_type, mc.cores = parallel::detectCores() - 3)
}

# Convert p-value results to a data frame for plotting
p_val_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  P_Value = unlist(p_val_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot p-values for each test type
# ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
#   geom_line() +  
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
#   labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests") +
#   theme_minimal() +
#   theme(legend.title = element_blank())

ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests using lognormal distribution under location shift") +
  theme_minimal() +
  theme(legend.title = element_blank())

# scale shift
# Define a function to perform the specified test for a given dimension
perform_test_lognorm_scale <- function(d, test_type) {
  # Generate X1 and X2 with the current dimension
  X1 <- MASS::mvrnorm(ss, rep(0, d), diag(1+sep, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
  X2 <- MASS::mvrnorm(ss, rep(0, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
  X1 <- exp(X1)
  X2 <- exp(X2)
  # Perform the specified test and return the p-value
  p_val <- test_fn(test_type, list(X1, X2), 0.05)
  return(p_val)
}

# Set the range of dimensions and test types
dimensions <- 2:200
test_types <- c("MCM", "MMCM", "FR")

# Initialize a list to store p-values for each test type
p_val_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  p_val_results[[test_type]] <- mclapply(dimensions, perform_test_lognorm_scale, test_type = test_type, mc.cores = parallel::detectCores() - 3)
}

# Convert p-value results to a data frame for plotting
p_val_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  P_Value = unlist(p_val_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot p-values for each test type
# ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
#   geom_line() +  
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
#   labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests") +
#   theme_minimal() +
#   theme(legend.title = element_blank())

ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests using lognormal distribution under scale shift") +
  theme_minimal() +
  theme(legend.title = element_blank())
#################################################
#################################################
# <<<t dist>>>
# location shift
# Define a function to perform the specified test for a given dimension
library(mvtnorm)
perform_test_t_loc <- function(d, test_type) {
  # Generate X1 and X2 with the current dimension
  X1 <- rmvt(ss, sigma = diag(1, d), df = 49, delta = rep(sep,d), type = "shifted")
  X2 <- rmvt(ss, sigma = diag(1, d), df = 49, delta = rep(0,d), type = "shifted")
  
  # Perform the specified test and return the p-value
  p_val <- test_fn(test_type, list(X1, X2), 0.05)
  return(p_val)
}

# Set the range of dimensions and test types
dimensions <- 2:200
test_types <- c("MCM", "MMCM", "FR")

# Initialize a list to store p-values for each test type
p_val_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  p_val_results[[test_type]] <- mclapply(dimensions, perform_test_t_loc, test_type = test_type, mc.cores = parallel::detectCores() - 3)
}

# Convert p-value results to a data frame for plotting
p_val_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  P_Value = unlist(p_val_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot p-values for each test type
# ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
#   geom_line() +  
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
#   labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests") +
#   theme_minimal() +
#   theme(legend.title = element_blank())

ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests using t distribution under location shift") +
  theme_minimal() +
  theme(legend.title = element_blank())

# scale shift
perform_test_t_scale <- function(d, test_type) {
  # Generate X1 and X2 with the current dimension
  X1 <- rmvt(ss, sigma = diag(1+sep, d), df = 49, delta = rep(0,d), type = "shifted")
  X2 <- rmvt(ss, sigma = diag(1, d), df = 49, delta = rep(0,d), type = "shifted")
  
  # Perform the specified test and return the p-value
  p_val <- test_fn(test_type, list(X1, X2), 0.05)
  return(p_val)
}

# Set the range of dimensions and test types
dimensions <- 2:200
test_types <- c("MCM", "MMCM", "FR")

# Initialize a list to store p-values for each test type
p_val_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  p_val_results[[test_type]] <- mclapply(dimensions, perform_test_t_scale, test_type = test_type, mc.cores = parallel::detectCores() - 3)
}

# Convert p-value results to a data frame for plotting
p_val_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  P_Value = unlist(p_val_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot p-values for each test type
# ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
#   geom_line() +  
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
#   labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests") +
#   theme_minimal() +
#   theme(legend.title = element_blank())

ggplot(p_val_df, aes(x = Log_Dimension, y = P_Value, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "P-Value", title = "P-Value vs Log(Dimension) for Different Tests using t distribution under scale shift") +
  theme_minimal() +
  theme(legend.title = element_blank())

#######################################################
#######################################################
# <<<Weibull dist>>>
# 