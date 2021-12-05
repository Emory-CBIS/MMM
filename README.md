---
output:
  pdf_document: default
  html_document: default
---
# Multi-Modal-Imaging

## About the Project
In this project, we develop a multimodal multilevel model (MMM) to infer brain maturation white matter structural connection and intrinsic functional connection from dMRI and resting-state fMRI. Scripts and data included here can be used to implement the method, reproduce the simulation studies and replicate the main findings in the manuscript. For the real data analysis, since we are not allowed to release the original data, we provide synthetic data generated based on the real data to replicate Figure 2-6 of the real data analysis in the paper. please note that the results based on the synthetic data may have some differences from the real data analysis results presented in our paper for the PNC study.

## Summary of Scripts in the Code Folder

### Model Related Scripts
+ model\_prep: script used to load library and all self-defined functions
+ MMM\_EM: function used to implement EM algorithm to obtain parameter estimation
+ calculate\_prob: function used to facilitate the probability calculation in the EM algorithm
+ simulate\_data: function used to simulate simulation data given certain parameter
+ initial\_param: function used to obtain initial estimation of the model 
+ initial\_random: function used to obtain initial estimation of the model (Random initialization for **$\alpha$** and **$\beta$** parameters)
+ utils: script with useful functions

### Simulation Related Scripts
+ simulation\_scripts: script to implement MMM for one simulation replicate
+ bootstrap\_scripts: script to implement MMM for one bootstrap replicate

### Real Data Analysis Related Scripts (included in RealDataAnalysis Folder)
+ analysis_prep: script used to load library and all self-defined functions for real data analysis
+ pnc_analysis_example: script used to implement MMM on 
+ pnc_analysis_boot: script used to run bootstrap and generate bootstrap replicates
+ pnc_analysis_results_summary: script used to summarize observed data and bootstrap data results to make statistical inference
+ analysis_help: script with useful functions to summarize results

## Procedure to Reproduce Simulation Data

1. Run model\_pre to load library and all functions needed
2. Run simulation\_scripts to generate one simulation replicate and save. Here are the detailed steps in simulation\_scripts.
   
   + Load fit.RData in the Data folder, which is used to help create parameters in the latent model in the simulation.
   + Simulate the Data
   + Set the initial parameter to estimation
   + Call EM_groupmod_log_Firth function in MMM_EM.r script to run MMM model
   + Save fitted model (fitsim.RData in the Data folder is an example of fitted model)
   
3. Run bootstrap\_scripts to generate multiple bootstrap replicates based on the results from Step 2 and then fit the model. Here are the detailed steps included in bootstrap\_scripts.

   + Load fitted model from Step 2 (for example, fitsim.RData) 
   + Obtain all the parameters needed to generate bootstrap data
   + Set the initial parameter to estimation
   + Call EM_groupmod_log_Firth function in MMM_EM.r script to run MMM model
   + Save fitted model for bootstrap data
   
4. Repeat Step 2 for 1000 times and for each simulation replicate and repeat Step 3 to obtain 200 bootstrap replicates.

## Procedure to Reproduce Real Data Analysis Results

1. Run analysis\_prep to load library and all functions needed
2. Run pnc\_analysis\_example to preprocess SC and FC data and fit MMM model
3. Run pnc\_analysis\_boot to generate multiple bootstrap replicates based on the results from 2 and then fit the model. 
4. Run pnc\_analysis\_results\_summary to analyze results and make statistical inference

## R version
R.4.1.1

## R pacakge used
Matrix_1.3-4; nor1mix_1.3-0; gtools_3.9.2; RcppArmadillo_0.10.7.0.0; Rcpp_1.0.7; MASS_7.3-54; nnet_7.3-16;
GGally_2.0.0; ggplot2_3.3.2; network_1.16.0; ComplexHeatmap_2.0.0; circlize_0.4.8; gplots_3.0.1.1;
