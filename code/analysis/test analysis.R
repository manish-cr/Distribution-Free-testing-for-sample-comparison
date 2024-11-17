setwd("/home/manish/Documents/NUS stuff/y4s1/DSA4288S-Mult/code/analysis")
library(tidyverse)
library(multicross)
library(FRmatch)
set.seed(1)
# there are 357 samples and 22278 genes. 2 classes
df <- read.csv("/home/manish/Documents/NUS stuff/y4s1/DSA4288S/data/liver_cancer.csv",
               stringsAsFactors = T)

# these will be the types
table(df['type'])
str(df[,1:3])

# train-test split
# set.seed(101)
# size <- floor(.8*nrow(df))
# sample_id <- sample.int(n=nrow(df),size=size)
# 
# train <- df[sample_id,]
# test <- df[-sample_id,]

# create a copy
df_copy <- data.frame(df)

# splitting into HCC and normal
HCC_indices <- which(df_copy$type=="HCC")
normal_indices <- which(df_copy$type == "normal")

# set data matrix
data_matrix_HCC <- as.matrix(df_copy[HCC_indices, -(1:2)])  # Exclude the first two columns (type and samples)
data_matrix_normal <- as.matrix(df_copy[normal_indices, -(1:2)])

# set the data list
data_list <- list(data_matrix_HCC,
                  data_matrix_normal)

p_val_mcm <- as.numeric(mcm(data_list, .05)[[1]])
p_val_mmcm <- as.numeric(mmcm(data_list, .05)[[1]])
p_val_fr <- FRtest(t(data_matrix_HCC), t(data_matrix_normal))[[5]]
