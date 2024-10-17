setwd("/home/manish/Documents/NUS stuff/y4s1/DSA4288S/code")
source("test_fn.R")
source("tree_gen.R")

# Define the GFS algorithm with dynamic tree creation and list handling
GFS <- function(features, alpha_thresh, data_list, test_fn) {
  
  # Create the tree dynamically
  tree <- create_tree(features)
  
  # Recursive function to perform the test and select the features
  recursive_test <- function(node, data_list, alpha_thresh, selected_features) {
    
    # Get the indices of the features in the current node
    current_features <- node$features
    
    # Create a sub-list for the current node's features
    data_sub <- lapply(data_list, function(data) data[, current_features, drop = FALSE])
    
    # Perform the test and get the p-value
    p_value <- test_fn(data_sub)
    
    # Adjusted p-value for the node
    adjusted_p_value <- p_value * (length(unlist(features)) / length(current_features))
    
    # If the adjusted p-value is greater than the threshold and the node has no children, stop
    if (adjusted_p_value > alpha_thresh && is.null(node$children)) {
      return(selected_features)  # No significant features
    }
    
    # If the node is a leaf and has a single feature
    if (adjusted_p_value <= alpha_thresh && length(current_features) == 1) {
      selected_features <- c(selected_features, current_features)
      return(selected_features)
    }
    
    # If the adjusted p-value is less than or equal to the threshold, check children nodes
    if (!is.null(node$children)) {
      for (child in node$children) {
        selected_features <- recursive_test(child, data_list, alpha_thresh, selected_features)
      }
    }
    
    return(selected_features)
  }
  
  # Start with the root node and an empty set of selected features
  selected_features <- recursive_test(tree, data_list, alpha_thresh, c())
  
  # Return the set of selected features
  return(selected_features)
}

# Example usage:
# Mock data matrices (replace with actual data)
X1 <- MASS::mvrnorm(100, rep(0, 4), diag(2, 4))
X1_signal <- rnorm(100,2,sd=2)
X1 <- cbind(X1_signal, X1)
X2 <- MASS::mvrnorm(100, rep(0, 5), diag(1, 5))
X3 <- MASS::mvrnorm(100, rep(0, 5), diag(3, 5))

# Define a data list
data_list <- list(X1, X2, X3)

# Use the MCM test function
test_fn <- function(data_sub) {
  return(as.numeric(mcm(data_sub, 0.05)[[1]]))
}

# Run the GFS algorithm with an alpha threshold of 0.05
alpha_thresh <- 0.05
features <- 1:5  # Example features
selected_features <- GFS(features, alpha_thresh, data_list, test_fn)

# Output the selected features
print(selected_features)
