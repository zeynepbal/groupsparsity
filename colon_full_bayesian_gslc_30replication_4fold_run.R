source("./classification/full_bayesian_group_sparse_linear_classification_variational_train.R")

library(AUC)
library(gglasso)

data(colon)
membership <- rep(1:20,each=5)

X_tra <- as.matrix(colon$x)
y_tra <- as.vector(colon$y)

D <- dim(X_tra)[2]

result_path <- "./results"
dir.create(result_path)

for (replication in 1:30) {
  negative_indices <- which(y_tra == -1)
  positive_indices <- which(y_tra == +1)
  
  #full Bayesian group sparse linear classification model parameters
  theta_alpha_list <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5)
  theta_beta <- 1
  lambda_alpha <- 1
  lambda_beta <- 1
  gamma_alpha <- 1
  gamma_beta <- 1
  eta_alpha <- 1
  eta_beta <- 1
  margin <- 1
  iteration <- 30
  progress <- 1
  
  fold_count <- 4 
  auroc_matrix_fbys <- matrix(NA, nrow = fold_count, ncol = length(theta_alpha_list))
  
  set.seed(12345 * replication)
  negative_allocation <- sample(rep(1:fold_count, ceiling(length(negative_indices) / fold_count)), length(negative_indices))
  positive_allocation <- sample(rep(1:fold_count, ceiling(length(positive_indices) / fold_count)), length(positive_indices))
  
  for (fold in 1:fold_count) {
    train_indices <- c(negative_indices[which(negative_allocation != fold)], positive_indices[which(positive_allocation != fold)])
    test_indices <- c(negative_indices[which(negative_allocation == fold)], positive_indices[which(positive_allocation == fold)])
    
    X_train <- X_tra[train_indices,]
    X_test <- X_tra[test_indices,]
    
    y_train <- y_tra[train_indices]
    y_test <- y_tra[test_indices]
    
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
      prediction_fbgslc <- fbgslc_model$bw$mu[1] + as.vector(X_test %*% (fbgslc_model$bw$mu[2:(D + 1)] * fbgslc_model$zeta$pi[membership])) 
      auroc_matrix_fbys[fold, t] <- auc(roc(prediction_fbgslc, as.factor(1 * (y_test == +1))))
      save("fbgslc_model", file = sprintf("%s/fb_repl_%d_fold_%d_prob_%g.RData", result_path, replication,fold, theta_alpha))
    }
    save("auroc_matrix_fbys", file = sprintf("%s/aurocfbys_repl_%d_prob_%g.RData", result_path, replication, theta_alpha))
    
  }
}

