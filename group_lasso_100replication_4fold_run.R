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
  lambda_grplasso <- c(seq(500, 1000, 100), 1100, 1200, 1300, 1400, 1500, 1900)
  
  fold_count <- 4 
  auroc_matrix_group_lasso <- matrix(NA, nrow = fold_count, ncol = length(lambda_grplasso)) 
  
  set.seed(12345 * replication)
  negative_allocation <- sample(rep(1:fold_count, ceiling(length(negative_indices) / fold_count)), length(negative_indices))
  positive_allocation <- sample(rep(1:fold_count, ceiling(length(positive_indices) / fold_count)), length(positive_indices))

  for (fold in 1:fold_count) {
    train_indices <- c(negative_indices[which(negative_allocation != fold)], positive_indices[which(positive_allocation != fold)])
    test_indices <- c(negative_indices[which(negative_allocation == fold)], positive_indices[which(positive_allocation == fold)])
    
    X_train <- X_tra[train_indices,]
    X_train_grp <- cbind(1, X_train)
    X_test <- X_tra[test_indices,]
    X_test_grp <- cbind(1, X_test)
    
    y_train <- y_tra[train_indices,]
    y_test <- y_tra[test_indices,]
    
    #Group Lasso model
    for (l in 1:length(lambda_grplasso)) {
      lambda <- lambda_grplasso[l]
      print(sprintf("running fold = %d, lambda = %g", fold, lambda))
      y_train_grp <- ifelse(y_train == -1, 0, y_train)
      y_test_grp <- ifelse(y_test == -1, 0, y_test)
      grplasso_model <- grplasso(x = X_train_grp, y = y_train_grp, model = LogReg(), lambda = lambda, index=membership_grp, center=FALSE, standardize = FALSE)
      prediction_grp <- predict(grplasso_model, newdata = as.data.frame(X_test_grp))
      auroc_matrix_group_lasso[fold, l] <- auc(roc(prediction_grp, as.factor(1 * (y_test_grp == +1))))
      save("grplasso_model", file = sprintf("%s/group_lasso_repl_%d_fold_%d_lambda_%g.RData", result_path, replication,fold, lambda))
    }
    save("auroc_matrix_group_lasso", file = sprintf("%s/auroc_group_lasso_repl_%d_lambda_%g.RData", result_path, replication, lambda))
  }
}