---
output:
  pdf_document: default
  html_document: default
---
# Multi-Modal-Imaging

## About the Project
In this project, we develop a multimodal multilevel model (MMM) to infer brain maturationin white matter structural connection and intrinsic functional connection from dMRI and resting-state fMRI. Scriptes and data attached here can be used to implement the model and replicate the simulation studies in the manuscript.

## Summary of Scripts in the Code Folder

+ model\_prep: function used to load library and all self-defined functions
+ MMM\_EM: function used to implement EM algorithm to obtain parameter estimation
+ calculate\_prob: function used to facilite the probabiliy calculation in the EM algorithm
+ simulate\_data: function used to simulate simulation data given cerain parameter
+ initial\_param: function used to obtain initial estimation of the model 
+ simluation\_scripts: script to implement MMM for one simulation replicate
+ bootstrap\_scripts: script to implement MMM for one bootstrap replicate

## Procedure to Reproduce Simulation Data

1. Run model\_pre to load library and all functions needed
2. Run simluation\_scripts to generate one simulation replicate and save. Here are the detailed steps in simluation\_scripts.
   
   + Load fit.RData in the Data folder, which is used to help create parameters in the latent model in the simuatlion.
   + Simulate the Data
   + Set the inital parameter to estimation
   + Call EM_groupmod_log_Firth function in MMM_EM.r script to run MMM model
   + Save fitted model (fit_1.RData in the Data folder is an example of fitted model)
   
3. Run bootstrap\_scripts to generate multiple bootstrap replicates based on the results from 2 and then fit the model. Here are the detalied steps in the Run bootstrap\_scripts.

   + Load fitted model from results 2 (for example fit_1.RData) 
   + Obtain all the parameters and Simulation bootstrap data
   + Set the inital parameter to estimation
   + Call EM_groupmod_log_Firth function in MMM_EM.r script to run MMM model
   + Save fitted model for boostrap data
   
4. Repeat 2 for 1000 times and for each simulation replicates repeat 3 for 200 times.

## R version
R.3.6.1

## R pacakge used
GGally_2.0.0; ggplot2_3.3.2; network_1.16.0; ComplexHeatmap_2.0.0; circlize_0.4.8;           
gplots_3.0.1.1; Matrix_1.2-18; nor1mix_1.3-0; gtools_3.8.1; RcppArmadillo_0.9.800.3.0;
Rcpp_1.0.3; MASS_7.3-51.5; nnet_7.3-12   