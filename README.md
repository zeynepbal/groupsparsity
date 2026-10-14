This repository contains R implementations of the Bayesian Group Sparse Linear Classification (BGSLC) and Full Bayesian Group Sparse Linear Classification (FBGSLC) algorithms presented in "Bringing group sparsity to Bayesian hierarchical models using multivariate Bernoulli distribution" by Zeynep Bal Türken, Gülay Başarır, and Mehmet Gönen. Both algorithms use a multivariate Bernoulli distribution for group variable selection within a Bayesian hierarchical framework and perform inference using variational approximation.

synthetic_data_demo.R file demonstrates how to train BGSLC and FBGSLC models in R using simulated data and compares the estimated coefficients with those of group Lasso.
bayesian_gslc_100replication_4fold_run.R file runs BGSLC on the splice-site detection dataset using 100 replications of stratified four-fold cross-validation.
full_bayesian_gslc_100replication_4fold_run.R file runs FBGSLC on the splice-site detection dataset using the same cross-validation procedure.
group_lasso_100replication_4fold_run.R file runs the group Lasso comparison on the splice-site detection dataset using the same cross-validation procedure.

BGSLC and FBGSLC methods
------------------------
* classification/bayesian_group_sparse_linear_classification_variational_train.R => variational training procedure for BGSLC in R
* classification/full_bayesian_group_sparse_linear_classification_variational_train.R => variational training procedure for FBGSLC in R

Splice-site detection data
--------------------------
* splice_site_detection/train_data_points.csv => nucleotide sequence predictors
* splice_site_detection/train_class_labels.csv => corresponding class labels

Requirements
------------
The code requires R and the arrangements, MASS, grplasso, and AUC packages. Install any missing packages in R using:

install.packages(c("arrangements", "MASS", "grplasso", "AUC"))

Running the code
----------------
Download or clone this repository and set the R working directory to the repository's main folder. Run the synthetic-data demonstration using:

source("synthetic_data_demo.R")

To reproduce the splice-site cross-validation experiments, run the following scripts separately from the same working directory:

source("bayesian_gslc_100replication_4fold_run.R")
source("full_bayesian_gslc_100replication_4fold_run.R")
source("group_lasso_100replication_4fold_run.R")
