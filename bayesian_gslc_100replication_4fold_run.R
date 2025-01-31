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

  #Bayesian group sparse linear classification model parameters
  lambda_alpha <- 1
  lambda_beta <- 1
  gamma_alpha <- 1
  gamma_beta <- 1
  p_list <-  c(seq(0.001, 0.005, 0.001), 0.01, 0.05, seq(0.1, 0.5, by = 0.1))
  margin <- 1
  iteration <- 100
  progress <- 1
  
  fold_count <- 4 
  auroc_matrix_bayesian <- matrix(NA, nrow = fold_count, ncol = length(p_list))
  
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
    
    #Bayesian group sparse linear classification model
    for (p in 1:length(p_list)) {
        prob <- p_list[p]
        print(sprintf("running fold = %d, p = %g", fold, prob))
        parameters <- list()
        parameters$lambda_alpha <- lambda_alpha
        parameters$lambda_beta <- lambda_beta
        parameters$gamma_alpha <- gamma_alpha
        parameters$gamma_beta <- gamma_beta
        parameters$p <- prob
        parameters$margin <- margin
        parameters$iteration <- iteration
        parameters$progress <- progress
        
        bgslc_model <- bayesian_group_sparse_linear_classification_variational_train(X = X_train, y = y_train, membership = membership, parameters = parameters)
        prediction_b <- bgslc_model$bw$mu[1] + as.vector(X_test %*% (bgslc_model$bw$mu[2:(D + 1)] * bgslc_model$zeta$pi[membership]))
        auroc_matrix_bayesian[fold, p] <- auc(roc(prediction_b, as.factor(1 * (y_test == +1))))
        save("bgslc_model", file = sprintf("%s/bayesian_repl_%d_fold_%d_prob_%g.RData", result_path, replication,fold, prob))
    }
    save("auroc_matrix_bayesian", file = sprintf("%s/auroc_bayesian_repl_%d_prob_%g.RData", result_path, replication, prob))
  }
}