# Function to create a hierarchical binary tree from features
create_tree <- function(features) {
  if (length(features) == 1) {
    return(list(features = features, children = NULL))  # Leaf node
  } else {
    mid <- ceiling(length(features) / 2)
    left_tree <- create_tree(features[1:mid])
    right_tree <- create_tree(features[(mid + 1):length(features)])
    return(list(
      features = features, 
      children = list(left_tree, right_tree)
    ))
  }
}

feat <- c(1:3)
tree <- create_tree(feat)
print(tree)
