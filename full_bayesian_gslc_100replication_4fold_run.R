source("./classification/bayesian_group_sparse_linear_classification_variational_train.R")
source("./classification/full_bayesian_group_sparse_linear_classification_variational_train.R")

library(AUC)
library(grplasso)

X_char <- read.csv("./splice_site_detection/train_data_points.csv", sep = ",", header = FALSE)
y_tra <- read.csv("./splice_site_detection/train_class_labels.csv", header = FALSE)

y_tra[y_tra == 2] <- 0
y_tra <- 2 * y_tra - 1

membership <- rep(c(1:7), each = 4)
G <- length(unique(membership))
N_train <- dim(X_char)[1]

X_tra <- matrix(0, N_train, G * 4)
for(j in 1:G) {
  X_tra[(X_char[,j] == "A"), (j - 1) * 4 + 1] <- 1
  X_tra[(X_char[,j] == "C"), (j - 1) * 4 + 2] <- 1
  X_tra[(X_char[,j] == "G"), (j - 1) * 4 + 3] <- 1
  X_tra[(X_char[,j] == "T"), (j - 1) * 4 + 4] <- 1
}

D <- dim(X_tra)[2]

result_path <- "./results"
dir.create(result_path)

for (replication in 1:100) {
  negative_indices <- which(y_tra == -1)
  positive_indices <- which(y_tra == +1)

  #full Bayesian group sparse linear classification model parameters
  theta_alpha_list <- c(10^seq(-5, -3, by = 1), seq(0.002, 0.009, 0.001), 0.01, 0.1, 1)
  theta_beta <- 1
  lambda_alpha <- 1
  lambda_beta <- 1
  gamma_alpha <- 1
  gamma_beta <- 1
  eta_alpha <- 1
  eta_beta <- 1
  margin <- 1
  iteration <- 100
  progress <- 1
  
  fold_count <- 4 
  auroc_matrix_full_bayesian <- matrix(NA, nrow = fold_count, ncol = length(theta_alpha_list))
  
  set.seed(12345 * replication)
  negative_allocation <- sample(rep(1:fold_count, ceiling(length(negative_indices) / fold_count)), length(negative_indices))
  positive_allocation <- sample(rep(1:fold_count, ceiling(length(positive_indices) / fold_count)), length(positive_indices))

  for (fold in 1:fold_count) {
    train_indices <- c(negative_indices[which(negative_allocation != fold)], positive_indices[which(positive_allocation != fold)])
    test_indices <- c(negative_indices[which(negative_allocation == fold)], positive_indices[which(positive_allocation == fold)])
    
    X_train <- X_tra[train_indices,]
    X_test <- X_tra[test_indices,]
    
    y_train <- y_tra[train_indices,]
    y_test <- y_tra[test_indices,]
    table(y_train)
    
    
    #full Bayesian group sparse linear classification model
    for (t in 1:length(theta_alpha_list)) {
      theta_alpha <- theta_alpha_list[t]
      print(sprintf("running fold = %d, theta = %g", fold, theta_alpha))
      parameters <- list()
      parameters$lambda_alpha <- lambda_alpha
      parameters$lambda_beta <- lambda_beta
      parameters$gamma_alpha <- gamma_alpha
      parameters$gamma_beta <- gamma_beta
      parameters$theta_alpha <- theta_alpha
      parameters$theta_beta <- theta_beta
      parameters$eta_alpha <- eta_alpha
      parameters$eta_beta <- eta_beta
      parameters$margin <- margin
      
      parameters$iteration <- iteration
      parameters$progress <- progress

      fbgslc_model <- full_bayesian_group_sparse_linear_classification_variational_train(X = X_train, y = y_train, membership = membership, parameters = parameters)
      prediction_f <- fbgslc_model$bw$mu[1] + as.vector(X_test %*% (fbgslc_model$bw$mu[2:(D + 1)] * fbgslc_model$zeta$pi[membership]))
      auroc_matrix_full_bayesian[fold, t] <- auc(roc(prediction_f, as.factor(1 * (y_test == +1))))
      save("fbgslc_model", file = sprintf("%s/full_bayesian_repl_%d_fold_%d_th_alpha_%g.RData", result_path, replication,fold, theta_alpha))
    }
    save("auroc_matrix_full_bayesian", file = sprintf("%s/auroc_full_bayesian_repl_%d_th_alpha_%g.RData", result_path, replication, theta_alpha))
  }
}