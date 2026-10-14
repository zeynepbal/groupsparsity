source("./classification/bayesian_group_sparse_linear_classification_variational_train.R")

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
  
  #Bayesian group sparse linear classification model parameters
  lambda_alpha <- 1
  lambda_beta <- 1
  gamma_alpha <- 1
  gamma_beta <- 1
  p_list <- c(0.001,0.005, seq(0.01, 0.05, by=0.01), seq(0.1, 0.5, by=0.1))
  margin <- 1
  iteration <- 100
  progress <- 1
  
  train_ratio_list <- 0.3
  
  auroc_b <- matrix(NA, 1, length(p_list))
  runtime_b <- matrix(NA, 1, length(p_list))
  
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
    
    #Bayesian group sparse linear classification model
    for (p in 1:length(p_list)) {
      prob <- p_list[p]
      print(sprintf("running replication = %d, train_ratio = %.2f, p = %g", replication, train_ratio, prob))
      parameters <- list()
      parameters$lambda_alpha <- 1
      parameters$lambda_beta <- 1
      parameters$gamma_alpha <- 1
      parameters$gamma_beta <- 1
      parameters$p <- prob
      parameters$margin <- margin
      parameters$iteration <- iteration
      parameters$progress <- progress
      
      t0 <- Sys.time()
      bgslc_model <- bayesian_group_sparse_linear_classification_variational_train(X = X_train, y = y_train, membership = membership, parameters = parameters)
      t1 <- Sys.time()
      
      prediction_b <- bgslc_model$bw$mu[1] + as.vector(X_test %*% (bgslc_model$bw$mu[2:(D + 1)] * bgslc_model$zeta$pi[membership]))
      
      auroc_b[p] <- auc(roc(prediction_b, as.factor(1 * (y_test == +1))))
      runtime_b[p] <- as.numeric(difftime(t1, t0, units = "secs"))
      
      save(bgslc_model, file = sprintf("%s/b_repl_%d_ratio_%g_prob_%g.RData", result_path, replication, train_ratio, prob))
    }
    save(auroc_b,runtime_b, file = sprintf("%s/bgslc_metrics_repl_%d_ratio_%g.RData", result_path, replication, train_ratio))
  }
}



