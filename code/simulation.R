set.seed(2)
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

A <- rbind(X1_dim_20,X2_dim_20)
source('test_gfs.R')
class_sizes <- c(100,100,100)
GFS(A,.05,class_sizes)