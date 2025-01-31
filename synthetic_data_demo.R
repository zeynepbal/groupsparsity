source("./classification/bayesian_group_sparse_linear_classification_variational_train.R")
source("./classification/full_bayesian_group_sparse_linear_classification_variational_train.R")
# create a synthetic data set
library(MASS)
set.seed(12345)
N <- 5000
coefs <- rep(0, 100)
membership <- c(rep(1:10, each = 10))
group_nonzeros <- c(10, 8, 6, 4, 2, 1)
for (g in 1:6){
  selected_features <- which(membership == g)
  coefs[sample(selected_features, size = group_nonzeros[g], replace = FALSE)] <- sample(c(-1, +1), size = group_nonzeros[g], replace = TRUE)
}
D <- length(coefs)

mu <- rep(0, D)
Sigma <- diag(0, D, D)
for (g in 1:max(membership)) {
  Sigma[which(membership == g), which(membership == g)] <- 0.25
}
diag(Sigma) <- 1.25
X <- mvrnorm(n = N, mu = mu, Sigma = Sigma)
y <- as.vector(X %*% coefs) + 4.0 * rnorm(N, mean = 0, sd = 1)
y <- ifelse(y > 0, +1, -1)

group_colors <- rep(c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6"), each = 10)


#ideal case
par(mar = c(6, 8, 0.9, 0.25), oma = c(0, 0, 0, 0))
barplot(coefs, las = 1, col = group_colors, xlab = "", cex.axis = 2)
mtext(side = 1, text = "Feature index", line = 3, cex = 2)
mtext(side = 2, text = "Feature weight", line = 6, cex = 2)

# Bayesian group sparse linear classification
parameters <- list()

parameters$lambda_alpha <- 1e-10
parameters$lambda_beta <- 1e+10
parameters$gamma_alpha <- 1
parameters$gamma_beta <- 1
parameters$margin <- 1
parameters$p <- 2 / 10

parameters$iteration <- 100
parameters$seed <- 123
parameters$progress <- 1

bgslc_model <- bayesian_group_sparse_linear_classification_variational_train(X = X, y = y, membership = membership, parameters = parameters)

par(mar = c(6, 8, 0.3, 0.25), oma = c(0, 0, 0, 0))
barplot(bgslc_model$bw$mu[-1] * bgslc_model$zeta$pi[membership], col = group_colors, 
        las = 1, xlab = "", cex.axis = 2)
mtext(side = 1, text = "Feature index", line = 3, cex = 2)
mtext(side = 2, text = "Feature weight", line = 6, cex = 2)


#Full Bayesian group sparse linear classification
parameters <- list()

parameters$lambda_alpha <- 1e-10
parameters$lambda_beta <- 1e+10
parameters$gamma_alpha <- 1
parameters$gamma_beta <- 1
parameters$theta_alpha <- 0.1
parameters$theta_beta <- 1
parameters$eta_alpha <- 1
parameters$eta_beta <- 1
parameters$margin <- 1


parameters$iteration <- 100
parameters$seed <- 123
parameters$progress <- 1

fbgslc_model <- full_bayesian_group_sparse_linear_classification_variational_train(X = X, y = y, membership = membership, parameters = parameters)

par(mar = c(6, 8, 0.3, 0.25), oma = c(0, 0, 0, 0))
barplot(fbgslc_model$bw$mu[-1] * fbgslc_model$zeta$pi[membership], col = group_colors, 
        las = 1, xlab = "", cex.axis = 2)
mtext(side = 1, text = "Feature index", line = 3, cex = 2)
mtext(side = 2, text = "Feature weight", line = 6, cex = 2)

#Group Lasso:
library(grplasso)
X <- cbind(1, X)
D <- dim(X)[2]

membership <- c(NA, rep(1:10, each = 10))

y_full_grp<- ifelse(y > 0, +1, 0)

grplasso_model <- grplasso(x = X, y = y_full_grp, model = LogReg(), lambda = 50, index=membership, center=FALSE, standardize = FALSE)

par(mar = c(6, 8, 0.3, 0.25), oma=c(0, 0, 0, 0))
barplot(grplasso_model$coefficients[,1], col = group_colors, 
        las = 1, xlab = "", cex.axis = 2)
mtext(side = 1, text = "Feature index", line = 3, cex = 2)
mtext(side = 2, text = "Feature weight", line = 6, cex = 2)

