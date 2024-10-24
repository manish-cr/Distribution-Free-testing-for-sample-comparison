set.seed(2)

# real signal set:
signal <- c(1:10)

# suppose that there are 100 observations
X1_signal <- MASS::mvrnorm(100,rep(1,10), diag(1,10))
X1_false_10_feat <- MASS::mvrnorm(100,rep(0,10), diag(1,10))
X1_false_20_feat <- MASS::mvrnorm(100,rep(0,20), diag(1,20))
X1_false_30_feat <- MASS::mvrnorm(100,rep(0,30), diag(1,30))


X2_signal <- MASS::mvrnorm(100,rep(5,10), diag(1,10))
X2_false_10_feat <- MASS::mvrnorm(100,rep(0,10), diag(1,10))
X2_false_20_feat <- MASS::mvrnorm(100,rep(0,20), diag(1,20))
X2_false_30_feat <- MASS::mvrnorm(100,rep(0,30), diag(1,30))

# doing the merge here
X1_dim_20 <- cbind(X1_signal, X1_false_10_feat)
X1_dim_30 <- cbind(X1_signal, X1_false_20_feat)
X1_dim_40 <- cbind(X1_signal, X1_false_30_feat)

X2_dim_20 <- cbind(X2_signal, X2_false_10_feat)
X2_dim_30 <- cbind(X2_signal, X2_false_20_feat)
X2_dim_40 <- cbind(X2_signal, X2_false_30_feat)

# simulation to find fwer and fdr
source('test_gfs.R')
data_list <- list(X1_dim_20, X2_dim_20)
features <- c(1:ncol(X1_dim_20))
alpha_thresh = .05
selected_features <- GFS(features, alpha_thresh, data_list, test_type = "MCM")

#==================================================
# Simulation approach (NEEDS CHECKING)
# Function to run GFS and compute FWER, FDR, and Power for each replicate
run_simulation <- function(n_dims, signal, alpha_thresh = 0.05, test_type = "MCM") {
  # Ensure randomness for each run
  set.seed(NULL)
  
  # Generate data for signal and false features
  X1_signal <- MASS::mvrnorm(100, rep(1, length(signal)), diag(1, length(signal)))
  X2_signal <- MASS::mvrnorm(100, rep(10, length(signal)), diag(1, length(signal)))
  
  false_features <- MASS::mvrnorm(100, rep(0, n_dims - length(signal)), diag(1, n_dims - length(signal)))
  X1 <- cbind(X1_signal, false_features)
  X2 <- cbind(X2_signal, false_features)
  
  # Prepare data for GFS
  data_list <- list(X1, X2)
  features <- c(1:ncol(X1))
  
  # Run GFS
  selected_features <- GFS(features, alpha_thresh, data_list, test_type = test_type)
  
  # Compute FWER, FDR, and Power
  S <- selected_features
  L <- signal
  
  # FWER: 1 if any false features (non-signal) are selected
  FWER <- as.numeric(any(setdiff(S, L)))
  
  # FDR: proportion of false discoveries (features in S that are not in L)
  FDR <- ifelse(length(S) > 0, length(setdiff(S, L)) / length(S), 0)
  
  # Power: 1 if all true features (L) are selected
  Power <- as.numeric(all(L %in% S))
  
  return(c(FWER = FWER, FDR = FDR, Power = Power))
}

# Function to simulate across replicates and compute mean FWER, FDR, and Power
simulate_over_reps <- function(n_dims, signal, n_reps = 100, alpha_thresh = 0.05, test_type = "MCM") {
  results <- replicate(n_reps, run_simulation(n_dims, signal, alpha_thresh, test_type))
  
  # Calculate the average for FWER, FDR, and Power
  mean_results <- rowMeans(results)
  return(mean_results)
}

# Example: Running the simulation for 20, 30, and 40 dimensions over 100 replicates
dims <- c(20, 30, 40)
dims <- c(50,60,70)
signal <- 1:10
n_reps <- 100

# Run the simulation for each dimension
sim_results <- lapply(dims, function(d) simulate_over_reps(d, signal, n_reps))
names(sim_results) <- paste0("Dim_", dims)

# Display results
sim_results

