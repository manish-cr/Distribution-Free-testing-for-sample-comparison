library(multicross)
library(FRmatch)
set.seed(1)
# Simulation Example when the user wants to test whether K=3 multivariate distributions are equal:
X1 = MASS::mvrnorm(100,rep(1,4),diag(2,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
X2 = MASS::mvrnorm(100,rep(0,4),diag(1,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
X3 = MASS::mvrnorm(10,rep(0,4),diag(3,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
mcm(list(X1,X2,X3),0.05)[[1]]
FRtest(t(X1),t(X2))

# using the painted turtles (baseline example)
# MCM test
# Creating vectors for males and females
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
# Assuming the males and females matrices are being used to represent two distributions
mcm(list(males, females), 0.05)

# MMCM test
mmcm(list(males, females), 0.05)

p_val <- FRtest(t(males), t(females))[[5]]
