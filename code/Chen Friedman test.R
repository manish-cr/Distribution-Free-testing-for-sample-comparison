# Load required library for distance calculation
library(igraph)

# Data: Males and Females matrices
males <- matrix(c(93, 74, 37,
                  94, 78, 35,
                  96, 80, 35,
                  101, 84, 39,
                  102, 85, 38,
                  103, 81, 37,
                  104, 83, 39,
                  106, 83, 39,
                  107, 82, 38,
                  112, 89, 40,
                  113, 88, 40,
                  114, 86, 40,
                  116, 90, 43,
                  117, 90, 41,
                  117, 91, 41,
                  119, 93, 41,
                  120, 89, 40,
                  120, 93, 44,
                  121, 95, 42,
                  125, 93, 45,
                  127, 96, 45,
                  128, 95, 45,
                  131, 95, 46,
                  135, 106, 47), ncol=3, byrow=TRUE)

females <- matrix(c(98, 81, 38,
                    103, 84, 38,
                    103, 86, 42,
                    105, 86, 40,
                    109, 88, 44,
                    123, 92, 50,
                    123, 95, 46,
                    133, 99, 51,
                    133, 102, 51,
                    133, 102, 51,
                    134, 100, 48,
                    136, 102, 49,
                    137, 98, 51,
                    138, 99, 51,
                    141, 105, 53,
                    147, 108, 57,
                    149, 107, 55,
                    153, 107, 56,
                    155, 115, 63,
                    155, 117, 60,
                    158, 115, 62,
                    159, 118, 63,
                    162, 124, 61,
                    177, 132, 67), ncol=3, byrow=TRUE)

# Another data test to check for non-rejection of H_0
X1 = MASS::mvrnorm(10,rep(0,4),diag(2,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
X2 = MASS::mvrnorm(10,rep(0,4),diag(1,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
data_list_1 <- list(males, females)
data_list_2 <- list(X1, X2)

# Function to calculate cosine similarity matrix
cosine_similarity <- function(X) {
  X_norm <- X / sqrt(rowSums(X^2))
  return(X_norm %*% t(X_norm))
}

# Function to calculate MST using cosine similarity
calculate_mst <- function(X) {
  sim_matrix <- cosine_similarity(X)
  # Convert similarity to distance for MST construction
  dist_matrix <- 1 - sim_matrix
  # Create a graph object from distance matrix, mode='lower' for undirected graph
  g <- graph_from_adjacency_matrix(dist_matrix)
  # Compute MST
  mst <- mst(g, weights = E(g)$weight)
  return(mst)
}

# Function to calculate R0, R1, R2
calculate_R_values <- function(mst, group_labels) {
  edges <- as_edgelist(mst)
  R0 <- 0; R1 <- 0; R2 <- 0
  
  for (edge in 1:nrow(edges)) {
    i <- edges[edge, 1]
    j <- edges[edge, 2]
    if (group_labels[i] != group_labels[j]) {
      R0 <- R0 + 1
    } else if (group_labels[i] == 0) {
      R1 <- R1 + 1
    } else {
      R2 <- R2 + 1
    }
  }
  return(list(R0=R0, R1=R1, R2=R2))
}

# Function to calculate the test-statistic with accordance to the formula
calculate_test_statistic <- function(R1, R2, mu1, mu2, cov_matrix, reg_lambda = 1e-5) {
  # Create the difference vector (R1 - mu1, R2 - mu2)
  diff_vector <- c(R1 - mu1, R2 - mu2)
  
  # Regularize the covariance matrix
  covariance_matrix_reg <- cov_matrix + diag(reg_lambda, nrow(cov_matrix))
  
  # Compute the inverse of the covariance matrix
  cov_inv <- solve(covariance_matrix_reg)
  
  # Compute the test statistic S
  S <- t(diff_vector) %*% cov_inv %*% diff_vector
  
  return(as.numeric(S))  # Return S as a scalar
}

# Calculate the MST on combined data
mst <- calculate_mst(combined_data)

# Get original R values
R_values <- calculate_R_values(mst, group_labels)
R1 <- R_values$R1
R2 <- R_values$R2

# Get the means (INCORRECT because R1 and R2 are just numbers here(must be done at permutation step))
mu1 <- mean(R1)
mu2 <- mean(R2)

# Estimate covariance matrix (this won't change during permutation)
cov_matrix <- calculate_covariance(R1_list = R1, R2_list = R2)

# Permutation testing function
chen_friedman <- function(data_list, num_permutations=1000) {
  # split up the data_list
  X1 <- data_list[[1]]
  X2 <- data_list[[2]]
  
  # combine
  combined_data <- rbind(X1, X2)
  group_labels <- c(rep(0, nrow(X1)), rep(1, nrow(X2)))
  
  original_mst <- calculate_mst(combined_data)
  original_R_values <- calculate_R_values(original_mst, group_labels)
  
  R1_values <- numeric(num_permutations)
  R2_values <- numeric(num_permutations)
  perm_stats <- numeric(num_permutations)
  
  # Perform permutations
  for (i in 1:num_permutations) {
    perm_labels <- sample(group_labels)
    perm_mst <- calculate_mst(combined_data)
    perm_R_values <- calculate_R_values(perm_mst, perm_labels)
    R1_values[i] <- perm_R_values$R1
    R2_values[i] <- perm_R_values$R2
  }
  
  # Calculate the means and covariance matrix
  mu1 <- mean(R1_values)
  mu2 <- mean(R2_values)
  cov_matrix <- cov(cbind(R1_values, R2_values))
  
  # Calculate test statistic for the original data
  original_statistic <- calculate_test_statistic(original_R_values$R1, original_R_values$R2, mu1, mu2, cov_matrix)
  
  # Calculate the test statistic for permutations
  for (i in 1:num_permutations) {
    perm_stats[i] <- calculate_test_statistic(R1_values[i], R2_values[i], mu1, mu2, cov_matrix)
  }
  
  # Calculate p-value
  p_value <- mean(perm_stats >= original_statistic)
  return(p_value)
}

# # Perform the permutation test
# p_value <- chen_friedman(data_list_1)
# 
# # Print the p-value
# print(paste("P-value:", p_value))
# 
# ############### non-rejection test
# # Combine X1 and X2
# 
# # Perform the permutation test
# p_value <- chen_friedman(data_list_2)
# 
# # Print the p-value
# print(paste("P-value:", p_value))