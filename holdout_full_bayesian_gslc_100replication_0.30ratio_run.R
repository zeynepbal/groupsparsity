source("./classification/full_bayesian_group_sparse_linear_classification_variational_train.R")

library(AUC)

X_char <- read.csv("train_data_points.csv", sep = ",", header = FALSE)
y_tra <- read.csv("train_class_labels.csv",header = FALSE)

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
  theta_alpha_list <- c(1e-04, 1e-03, 1e-02, seq(0.02, 0.09, by = 0.01), 1)
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
  
  train_ratio_list <- 0.3

  auroc_fb <- matrix(NA, 1, length(theta_alpha_list))
  runtime_fb <- matrix(NA, 1, length(theta_alpha_list))
  
  for (train_ratio  in train_ratio_list) {
    set.seed(12345 * replication)
    train_neg <- sample(negative_indices, floor(length(negative_indices) * train_ratio))
    train_pos <- sample(positive_indices, floor(length(positive_indices) * train_ratio))
    
    test_neg <- setdiff(negative_indices, train_neg)
    test_pos <- setdiff(positive_indices, train_pos)
    
    train_indices <- c(train_neg, train_pos)
    test_indices <- c(test_neg, test_pos)
    
    X_train <- X_tra[train_indices,]
    X_train_grp <- cbind(1, X_train)
    X_test <- X_tra[test_indices,]
    X_test_grp <- cbind(1, X_test)
    
    y_train <- y_tra[train_indices,]
    y_test <- y_tra[test_indices,]
    
    #full Bayesian group sparse linear classification model
    for (t in 1:length(theta_alpha_list)) {
      theta_alpha <- theta_alpha_list[t]
      print(sprintf("running replication = %d, train_ratio = %.2f, theta_alpha = %g", replication, train_ratio, theta_alpha))
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
      
      t0 <- Sys.time()
      fbgslc_model <- full_bayesian_group_sparse_linear_classification_variational_train(X = X_train, y = y_train, membership = membership, parameters = parameters)
      t1 <- Sys.time()
      
      prediction_f <- fbgslc_model$bw$mu[1] + as.vector(X_test %*% (fbgslc_model$bw$mu[2:(D + 1)] * fbgslc_model$zeta$pi[membership]))
      
      auroc_fb[t] <- auc(roc(prediction_f, as.factor(1 * (y_test == +1))))
      runtime_fb[t] <- as.numeric(difftime(t1, t0, units = "secs"))
      
      save(fbgslc_model, file = sprintf("%s/f_repl_%d_ratio_%.2f_th_alpha_%g.RData", result_path, replication,train_ratio, theta_alpha))
    }
    save(auroc_fb,runtime_fb, file = sprintf("%s/fbgslc_metrics_repl_%d_ratio_%.2f.RData", result_path, replication, train_ratio))
  }
}



