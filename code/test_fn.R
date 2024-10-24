library(multicross)
library(FRmatch)
setwd("/home/manish/Documents/NUS stuff/y4s1/DSA4288S-Mult/code")
source("Chen Friedman test.R")

# Updated test function that directly takes a data list (list of matrices)
test_fn <- function(test_type, data_list, level = 0.05) {
  if (test_type == "MCM") {
    p_val <- as.numeric(mcm(data_list, level)[[1]])
    return(p_val)
  } else if (test_type == "MMCM") {
    p_val <- as.numeric(mmcm(data_list, level)[[1]])
    return(p_val)
  } else if (test_type == "CF"){
    p_val <- as.numeric(chen_friedman(data_list))
    return(p_val)
  }
}
