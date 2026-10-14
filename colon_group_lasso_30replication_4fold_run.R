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
  
  #Group Lasso model parameters:
  lambda_list <-  c(0.0001,0.001, 0.005, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.03, 0.04, 0.05)
  fold_count <- 4 
  
  auroc_matrix_gg <- matrix(NA, nrow = fold_count, ncol = length(lambda_list))
  
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
    
    #Group Lasso model
    for (l in 1:length(lambda_list)) {
      lmd <- lambda_list[l]
      print(sprintf("running fold = %d, p = %g", fold, lmd))
      gglasso_model <- gglasso(x=X_train, y=y_train, group=membership, loss="logit", lambda=lmd)
      prediction_g <- predict(gglasso_model, newx = X_test)
      auroc_matrix_gg[fold, l] <- auc(roc(prediction_g, as.factor(1 * (y_test == +1))))
      save("gglasso_model", file = sprintf("%s/g_repl_%d_fold_%d_prob_%g.RData", result_path, replication,fold, lmd))
    }
    save("auroc_matrix_gg", file = sprintf("%s/aurocgg_repl_%d_prob_%g.RData", result_path, replication, lmd))
    
  }
}
