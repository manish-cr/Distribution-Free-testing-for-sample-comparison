setwd("/home/manish/Documents/NUS stuff/y4s1/DSA4288S-Mult/code/analysis")
source("../test_fn.R")
library(tidyverse)
library(multicross)
library(FRmatch)

# Sample size and list of dimensions
ss <- 50
dimensions <- 2:200
test_types <- c("MCM", "MMCM", "FR")
num_simulations <- 20  # Number of repeated simulations for each dimension

# <<<Normal distribution>>>
# Function to generate samples and calculate power based on repeated tests
calculate_power_normal_loc <- function(d, sep, test_type) {
  # Run simulations in parallel and count rejections
  rejections <- parallel::mclapply(1:num_simulations, function(i) {
    # Generate X1 and X2 with location shift for the specified dimension
    X1 <- MASS::mvrnorm(ss, rep(sep, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
    X2 <- MASS::mvrnorm(ss, rep(0, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
    
    # Perform the specified test and check if p-value < 0.05
    p_val <- test_fn(test_type, list(X1, X2), 0.05)
    return(as.numeric(p_val < 0.05))  # Return 1 if reject, 0 otherwise
  }, mc.cores = parallel::detectCores() - 3)
  
  # Calculate power as proportion of rejections
  power <- mean(unlist(rejections))
  return(power)
}

# Initialize a list to store power results for each test type
power_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  power_results[[test_type]] <- mclapply(dimensions, function(d) calculate_power_normal_loc(d, sep = 0.5, test_type = test_type),
                                         mc.cores = parallel::detectCores() - 3)
}

# Convert power results to a data frame for plotting
power_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  Power = unlist(power_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot power vs log dimension for each test type
ggplot(power_df, aes(x = Log_Dimension, y = Power, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "Power", title = "Power vs Log(Dimension) for Different Tests using normal distribution with Location Shift") +
  theme_minimal() +
  theme(legend.title = element_blank())
##########################################
##########################################
# scale shift
calculate_power_normal_scale <- function(d, sep, test_type) {
  # Run simulations in parallel and count rejections
  rejections <- parallel::mclapply(1:num_simulations, function(i) {
    # Generate X1 and X2 with location shift for the specified dimension
    X1 <- MASS::mvrnorm(ss, rep(0, d), diag(1+sep, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
    X2 <- MASS::mvrnorm(ss, rep(0, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
    
    # Perform the specified test and check if p-value < 0.05
    p_val <- test_fn(test_type, list(X1, X2), 0.05)
    return(as.numeric(p_val < 0.05))  # Return 1 if reject, 0 otherwise
  }, mc.cores = parallel::detectCores() - 3)
  
  # Calculate power as proportion of rejections
  power <- mean(unlist(rejections))
  return(power)
}

# Initialize a list to store power results for each test type
power_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  power_results[[test_type]] <- mclapply(dimensions, function(d) calculate_power_normal_scale(d, sep = 0.5, test_type = test_type),
                                         mc.cores = parallel::detectCores() - 3)
}

# Convert power results to a data frame for plotting
power_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  Power = unlist(power_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot power vs log dimension for each test type
ggplot(power_df, aes(x = Log_Dimension, y = Power, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "Power", title = "Power vs Log(Dimension) for Different Tests using normal distribution with Scale Shift") +
  theme_minimal() +
  theme(legend.title = element_blank())
######################################
######################################
# <<<lognormal distribution>>>
calculate_power_lognormal_loc <- function(d, sep, test_type) {
  # Run simulations in parallel and count rejections
  rejections <- parallel::mclapply(1:num_simulations, function(i) {
    # Generate X1 and X2 with location shift for the specified dimension
    X1 <- MASS::mvrnorm(ss, rep(sep, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
    X2 <- MASS::mvrnorm(ss, rep(0, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
    X1 <- exp(X1)
    X2 <- exp(X2)
    
    # Perform the specified test and check if p-value < 0.05
    p_val <- test_fn(test_type, list(X1, X2), 0.05)
    return(as.numeric(p_val < 0.05))  # Return 1 if reject, 0 otherwise
  }, mc.cores = parallel::detectCores() - 3)
  
  # Calculate power as proportion of rejections
  power <- mean(unlist(rejections))
  return(power)
}

# Initialize a list to store power results for each test type
power_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  power_results[[test_type]] <- mclapply(dimensions, function(d) calculate_power_lognormal_loc(d, sep = 0.5, test_type = test_type),
                                         mc.cores = parallel::detectCores() - 3)
}

# Convert power results to a data frame for plotting
power_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  Power = unlist(power_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot power vs log dimension for each test type
ggplot(power_df, aes(x = Log_Dimension, y = Power, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "Power", title = "Power vs Log(Dimension) for Different Tests using lognormal distribution with Location Shift") +
  theme_minimal() +
  theme(legend.title = element_blank())
#######################################
#######################################
# Scale shift
calculate_power_lognormal_scale <- function(d, sep, test_type) {
  # Run simulations in parallel and count rejections
  rejections <- parallel::mclapply(1:num_simulations, function(i) {
    # Generate X1 and X2 with location shift for the specified dimension
    X1 <- MASS::mvrnorm(ss, rep(0, d), diag(1+sep, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
    X2 <- MASS::mvrnorm(ss, rep(0, d), diag(1, d), tol=1e-6, empirical=FALSE, EISPACK=FALSE)
    X1 <- exp(X1)
    X2 <- exp(X2)
    
    # Perform the specified test and check if p-value < 0.05
    p_val <- test_fn(test_type, list(X1, X2), 0.05)
    return(as.numeric(p_val < 0.05))  # Return 1 if reject, 0 otherwise
  }, mc.cores = parallel::detectCores() - 3)
  
  # Calculate power as proportion of rejections
  power <- mean(unlist(rejections))
  return(power)
}

# Initialize a list to store power results for each test type
power_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  power_results[[test_type]] <- mclapply(dimensions, function(d) calculate_power_lognormal_scale(d, sep = 0.5, test_type = test_type),
                                         mc.cores = parallel::detectCores() - 3)
}

# Convert power results to a data frame for plotting
power_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  Power = unlist(power_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot power vs log dimension for each test type
ggplot(power_df, aes(x = Log_Dimension, y = Power, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "Power", title = "Power vs Log(Dimension) for Different Tests using lognormal distribution with Scale Shift") +
  theme_minimal() +
  theme(legend.title = element_blank())
############################
############################
# <<<t dist>>>
library(mvtnorm)
calculate_power_t_loc <- function(d, sep, test_type) {
  # Run simulations in parallel and count rejections
  rejections <- parallel::mclapply(1:num_simulations, function(i) {
    # Generate X1 and X2 with location shift for the specified dimension
    X1 <- rmvt(ss, sigma = diag(1, d), df = 49, delta = rep(sep,d), type = "shifted")
    X2 <- rmvt(ss, sigma = diag(1, d), df = 49, delta = rep(0,d), type = "shifted")
    
    # Perform the specified test and check if p-value < 0.05
    p_val <- test_fn(test_type, list(X1, X2), 0.05)
    return(as.numeric(p_val < 0.05))  # Return 1 if reject, 0 otherwise
  }, mc.cores = parallel::detectCores() - 3)
  
  # Calculate power as proportion of rejections
  power <- mean(unlist(rejections))
  return(power)
}

# Initialize a list to store power results for each test type
power_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  power_results[[test_type]] <- mclapply(dimensions, function(d) calculate_power_t_loc(d, sep = 0.5, test_type = test_type),
                                         mc.cores = parallel::detectCores() - 3)
}

# Convert power results to a data frame for plotting
power_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  Power = unlist(power_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot power vs log dimension for each test type
ggplot(power_df, aes(x = Log_Dimension, y = Power, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "Power", title = "Power vs Log(Dimension) for Different Tests using t distribution with Location Shift") +
  theme_minimal() +
  theme(legend.title = element_blank())
#######################################
#######################################
calculate_power_t_scale <- function(d, sep, test_type) {
  # Run simulations in parallel and count rejections
  rejections <- parallel::mclapply(1:num_simulations, function(i) {
    # Generate X1 and X2 with location shift for the specified dimension
    X1 <- rmvt(ss, sigma = diag(1+sep, d), df = 49, delta = rep(0,d), type = "shifted")
    X2 <- rmvt(ss, sigma = diag(1, d), df = 49, delta = rep(0,d), type = "shifted")
    
    # Perform the specified test and check if p-value < 0.05
    p_val <- test_fn(test_type, list(X1, X2), 0.05)
    return(as.numeric(p_val < 0.05))  # Return 1 if reject, 0 otherwise
  }, mc.cores = parallel::detectCores() - 3)
  
  # Calculate power as proportion of rejections
  power <- mean(unlist(rejections))
  return(power)
}

# Initialize a list to store power results for each test type
power_results <- list()

# Run tests in parallel for each test type and each dimension
for (test_type in test_types) {
  power_results[[test_type]] <- mclapply(dimensions, function(d) calculate_power_t_loc(d, sep = 0.5, test_type = test_type),
                                         mc.cores = parallel::detectCores() - 3)
}

# Convert power results to a data frame for plotting
power_df <- data.frame(
  Dimension = rep(dimensions, times = length(test_types)),
  Log_Dimension = rep(log(dimensions), times = length(test_types)),
  Power = unlist(power_results),
  Test_Type = rep(test_types, each = length(dimensions))
)

# Plot power vs log dimension for each test type
ggplot(power_df, aes(x = Log_Dimension, y = Power, color = Test_Type)) +
  geom_line(alpha = 0.3) +  # Make the original line semi-transparent for clarity
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +  # Add a smoothed line with LOESS
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # Y-axis ticks at intervals of 0.1
  labs(x = "Log(Dimension)", y = "Power", title = "Power vs Log(Dimension) for Different Tests using t distribution with Location Shift") +
  theme_minimal() +
  theme(legend.title = element_blank())