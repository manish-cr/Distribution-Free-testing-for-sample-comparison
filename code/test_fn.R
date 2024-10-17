library(multicross)

# Updated test function that directly takes a data list (list of matrices)
test_fn <- function(test_type = "MCM", data_list, level = 0.05) {
  if (test_type == "MCM") {
    p_val <- as.numeric(mcm(data_list, level)[[1]])
    return(p_val)
  } else if (test_type == "MMCM") {
    p_val <- as.numeric(mmcm(data_list, level)[[1]])
    return(p_val)
  }
}
