This repository contains R implementations of the Bayesian Group Sparse Linear Classification (BGSLC) and Full Bayesian Group Sparse Linear Classification (FBGSLC) algorithms presented in "Bringing group sparsity to Bayesian hierarchical models using multivariate Bernoulli distribution" by Zeynep Bal Türken, Gülay Başarır, and Mehmet Gönen. Both algorithms use a multivariate Bernoulli distribution for group variable selection within a Bayesian hierarchical framework and perform inference using variational approximation.

synthetic_data_demo.R file demonstrates how to train BGSLC and FBGSLC models in R using simulated data and compares the estimated coefficients with those of group lasso.  
bayesian_gslc_100replication_4fold_run.R file runs BGSLC on the splice-site detection dataset using 100 replications of stratified four-fold cross-validation.  
full_bayesian_gslc_100replication_4fold_run.R file runs FBGSLC on the splice-site detection dataset using the same cross-validation procedure.  
group_lasso_100replication_4fold_run.R file runs the group Lasso comparison on the splice-site detection dataset using the same cross-validation procedure.  
holdout_bayesian_gslc_100replication_0.30ratio_run.R file evaluates BGSLC on the splice-site detection dataset using 100 stratified holdout replications with 30% training and 70% test data.  
holdout_full_bayesian_gslc_100replication_0.30ratio_run.R file evaluates FBGSLC using the same holdout procedure.  
holdout_group_lasso_100replication_0.30ratio_run.R file evaluates group Lasso using the same holdout procedure.  
colon_bayesian_gslc_30replication_4fold_run.R file evaluates BGSLC on the colon gene expression dataset using 30 replications of stratified four-fold cross-validation.  
colon_full_bayesian_gslc_30replication_4fold_run.R file evaluates FBGSLC on the colon dataset using the same cross-validation procedure.  
colon_group_lasso_30replication_4fold_run.R file evaluates group Lasso on the colon dataset using the same cross-validation procedure.

BGSLC and FBGSLC methods
------------------------
* classification/bayesian_group_sparse_linear_classification_variational_train.R => variational training procedure for BGSLC in R
* classification/full_bayesian_group_sparse_linear_classification_variational_train.R => variational training procedure for FBGSLC in R

Splice-site detection data
--------------------------
* splice_site_detection/train_data_points.csv => nucleotide sequence predictors
* splice_site_detection/train_class_labels.csv => corresponding class labels

Colon gene expression data
--------------------------
The colon dataset is available in the gglasso R package and is loaded by the colon experiment scripts using `data(colon)`.

Requirements
------------
The code requires R and the arrangements, MASS, grplasso, gglasso, and AUC packages. Install any missing packages in R using:

```r
install.packages(c("arrangements", "MASS", "grplasso", "gglasso", "AUC"))
```

Running the code
----------------
Download or clone this repository and set the R working directory to the repository's main folder. Run the synthetic-data demonstration using:

```r
source("synthetic_data_demo.R")
```

To reproduce the splice-site cross-validation experiments, run the following scripts separately from the same working directory:

```r
source("bayesian_gslc_100replication_4fold_run.R")
source("full_bayesian_gslc_100replication_4fold_run.R")
source("group_lasso_100replication_4fold_run.R")
```

To run the splice-site holdout experiments (100 replications, 30% training and 70% test data), use:

```r
source("holdout_bayesian_gslc_100replication_0.30ratio_run.R")
source("holdout_full_bayesian_gslc_100replication_0.30ratio_run.R")
source("holdout_group_lasso_100replication_0.30ratio_run.R")
```

To run the colon experiments (30 replications of stratified four-fold cross-validation), use:

```r
source("colon_bayesian_gslc_30replication_4fold_run.R")
source("colon_full_bayesian_gslc_30replication_4fold_run.R")
source("colon_group_lasso_30replication_4fold_run.R")
```
