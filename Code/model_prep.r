# source package
library(nnet)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(gtools)
library(nor1mix)
library(Matrix)

# source useful function
# set path directory
path <- paste(getwd(), "/Code/", sep = "")
source(paste(path, "MMM_EM.r", sep = ""))
source(paste(path, "brmultilogit.r", sep = ""))
source(paste(path, "calculate_prob.r", sep = ""))
source(paste(path, "initial_param.r", sep = ""))
source(paste(path, "initial_random.r", sep = ""))
source(paste(path, "simulate_data.r", sep = ""))
source(paste(path, "utils.r", sep = ""))
sourceCpp(paste(path, "fastEStep.cpp", sep = "")) 


