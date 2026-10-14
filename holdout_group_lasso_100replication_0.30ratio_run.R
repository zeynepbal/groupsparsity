library(AUC)
library(grplasso)

X_char <- read.csv("train_data_points.csv", sep = ",", header = FALSE)
y_tra <- read.csv("train_class_labels.csv",header = FALSE)

y_tra[y_tra == 2] <- 0
y_tra <- 2 * y_tra - 1

membership <- rep(c(1:7), each = 4)
G <- length(unique(membership))
N_train <- dim(X_char)[1]

membership_grp <- c(NA, rep(c(1:7), each = 4)) #for Group Lasso

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
  
  #Group Lasso model parameters:
  lambda_grplasso <- c(200, 250, 300, 350, 400, 450, 500, 550, seq(600, 1000, 100))
  
  train_ratio_list <- 0.3
  
  auroc_g <- matrix(NA, 1, length(lambda_grplasso))
  runtime_g <- matrix(NA, 1, length(lambda_grplasso))
  
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
    
    #Group Lasso model
    for (l in 1:length(lambda_grplasso)) {
      lambda <- lambda_grplasso[l]
      print(sprintf("running replication = %d, train_ratio = %.2f, lambda = %g", replication, train_ratio, lambda))
      y_train_grp <- ifelse(y_train==-1,0, y_train)
      y_test_grp <- ifelse(y_test==-1,0, y_test)
      
      t0 <- Sys.time()
      grplasso_model <- grplasso(x= X_train_grp, y= y_train_grp, model = LogReg(), lambda = lambda, index=membership_grp, center=FALSE, standardize = FALSE)
      t1 <- Sys.time()
      
      prediction_grp <- predict(grplasso_model, newdata = as.data.frame(X_test_grp))
      
      auroc_g[l] <- auc(roc(prediction_grp, as.factor(1 * (y_test_grp == +1))))
      runtime_g[l] <- as.numeric(difftime(t1, t0, units = "secs"))
      
      save(grplasso_model, file = sprintf("%s/g_repl_%d_ratio_%.1f_lambda_%g.RData", result_path, replication,train_ratio, lambda))
    }
    save(auroc_g, runtime_g, file = sprintf("%s/grp_metrics_repl_%d_ratio_%g.RData", result_path, replication, train_ratio))
  }
}



