library(multicross)
library(FRmatch)

# Updated test function that directly takes a data list (list of matrices)
test_fn <- function(test_type, data_list, level = 0.05) {
  if (test_type == "MCM") {
    p_val <- as.numeric(mcm(data_list, level)[[1]])
    return(p_val)
  } else if (test_type == "MMCM") {
    p_val <- as.numeric(mmcm(data_list, level)[[1]])
    return(p_val)
  } else if(test_type == "FR"){
    p_val <- as.numeric(FRtest(data_list[[1]],data_list[[2]])[[5]])
    if(is.na(p_val)){ # this is used because p_val tends to be NAN for small values of p_val
      p_val <- .01
    }
    return(p_val)
  }
}
