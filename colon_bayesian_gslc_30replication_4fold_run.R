source("classification/bayesian_group_sparse_linear_classification_variational_train.R")

library(AUC)
library(grplasso)
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
  
  #Bayesian group sparse linear classification model parameters
  lambda_alpha <- 1
  lambda_beta <- 1
  gamma_alpha <- 1
  gamma_beta <- 1
  p_list <- seq(0.1, 0.9, by = 0.1)
  margin <- 1
  iteration <- 30
  progress <- 1
  
  fold_count <- 4 
  auroc_matrix_bys <- matrix(NA, nrow = fold_count, ncol = length(p_list))
  
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
    
    #Bayesian group sparse linear classification model
    for (p in 1:length(p_list)) {
      prob <- p_list[p]
      print(sprintf("running fold = %d, p = %g", fold, prob))
      parameters <- list()
      parameters$lambda_alpha <- 1
      parameters$lambda_beta <- 1
      parameters$gamma_alpha <- 1
      parameters$gamma_beta <- 1
      parameters$p <- prob
      parameters$margin <- margin
      
      parameters$iteration <- iteration
      parameters$progress <- progress
      bgslc_model <- bayesian_group_sparse_linear_classification_variational_train(X = X_train, y = y_train, membership = membership, parameters = parameters)
      prediction_b <- bgslc_model$bw$mu[1] + as.vector(X_test %*% (bgslc_model$bw$mu[2:(D + 1)] * bgslc_model$zeta$pi[membership]))
      auroc_matrix_bys[fold, p] <- auc(roc(prediction_b, as.factor(1 * (y_test == +1))))
      save("bgslc_model", file = sprintf("%s/b_repl_%d_fold_%d_prob_%g.RData", result_path, replication,fold, prob))
    }
    save("auroc_matrix_bys", file = sprintf("%s/aurocbys_repl_%d_prob_%g.RData", result_path, replication, prob))
  }
}

