setwd("/home/manish/Documents/NUS stuff/y4s1/DSA4288S-Mult/code/analysis")
source("../test_fn.R")
library(MASS)
library(parallel)
library(dplyr)
library(knitr)
library(kableExtra)
set.seed(2)
# Sample size and list of dimensions, scale shifts, and test types
ss <- 50
dimensions <- c(5, 10, 50, 100, 200, 300, 500)
scale_shifts <- c(0.4, 0.9, 1.4, 1.9, 2.4)
test_types <- c("MCM", "MMCM", "FR")
num_simulations <- 50  # Number of simulations to calculate power

# Function to perform the lognormal test and calculate power for given parameters
calculate_power <- function(d, sep, test_type) {
  # Run simulations in parallel and count rejections
  rejections <- parallel::mclapply(1:num_simulations, function(i) {
    # Generate X1 and X2 with the specified scale shift and dimension
    X1 <- MASS::mvrnorm(ss, rep(0, d), diag(1 + sep, d), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    X2 <- MASS::mvrnorm(ss, rep(0, d), diag(1, d), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    X1 <- exp(X1)
    X2 <- exp(X2)
    
    # Perform the specified test and check if p-value < 0.05
    p_val <- test_fn(test_type, list(X1, X2), 0.05)
    return(as.numeric(p_val < 0.05))  # Return 1 if reject, 0 otherwise
  }, mc.cores = parallel::detectCores() - 3)
  
  # Calculate power as the proportion of rejections
  power <- mean(unlist(rejections))
  return(power)
}

# Initialize list to store power results for each combination of scale shift, dimension, and test type
power_results <- list()

# Parallelize across scale shifts and test types
for (sep in scale_shifts) {
  for (test_type in test_types) {
    # For each combination of scale shift and test type, compute power across dimensions
    power_vals <- mclapply(dimensions, function(d) calculate_power(d, sep, test_type),
                           mc.cores = parallel::detectCores() - 3)
    # Store results in the list
    power_results[[paste(test_type, sep)]] <- unlist(power_vals)
  }
}

# Convert the results to a data frame for tabular display
power_df <- data.frame(
  Dimension = rep(dimensions, times = length(scale_shifts) * length(test_types)),
  Scale_Shift = rep(scale_shifts, each = length(dimensions) * length(test_types)),
  Test_Type = rep(test_types, each = length(dimensions), times = length(scale_shifts)),
  Power = unlist(power_results)
)

# Reshape data to a wide-format table for a clean display
library(reshape2)
table_df <- dcast(power_df, Scale_Shift + Test_Type ~ Dimension, value.var = "Power")

# Print the table for viewing
print(table_df)

# Create a multi-index table with `kable` for better presentation
kable(table_df, row.names = FALSE, align = "c", caption = "Power for Different Scale Shifts and Dimensions") %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  add_header_above(c(" " = 2, "Dimensions" = length(dimensions)))
