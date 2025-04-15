# HDFLM-fEnet

This repo contains .R funcs and .rds data for simulating high-dimensional functional data and estimating the functional coefficients using our proposed functional elastic-net (fEnet) algorithms. The preprint of the paper can be find on ArXiv at https://arxiv.org/pdf/2310.14419. Simulation studies for For Scenario I, (n, p, q) = (100, 200, 10), rho = 0.75 can be reproduced using the R code. 

## Overview

- **simu_funcs.R**: This file contains functions to simulate high-dimensional functional data.

- **algo_enet_representer.R**: This file contains functions to estimate the functional coefficients using fEnet. Note that the provided code are based on a specific RKHS describe in the simulation of the paper. The user may modify the *K_sqrt* function (i.e. the squared root of the reproduce kernel) as needed.  

- **algo_refine.R**: This file contains the function to compute refine estimator the functional coefficients using ridge regression.

- **algo_scad.R**: This file contains function to estimate the functional coefficients using FLM-SCAD (Xue and Yao, 2021).

- **example_fenet.R**: This file contains an example to simulate data and fit the data using fenet and refined estimator (run this file to check the algorithm).

- **example_scad.R**: This file contains an example to simulate data and fit the data using FLM-SCAD.

- **cv_enet2_auto.R**: This file is used to perform 200 replicates of data simulation, estimation, and prediction using fEnet (with grid search of the optimal hyperparameters). This file may take 1-2 hours for 10 CPU cores.

- **cv_scad2_auto.R**: This file is used to perform 200 replicates of data simulation, estimation, and prediction using FLM-SCAD (with grid search of the optimal hyperparameters). This file may take 1-2 hours for 10 CPU cores.

- **cv_enet2_auto.rds**: This RDS file contains results for FPR, FNR, RMSPE, MND, RER of the 200 replicates of estimation and prediction results (fEnet) for all combinations of hyperparameters.

- **cv_scad2_auto.rds**: This RDS file contains results for FPR, FNR, RMSPE, MND, RER of the 200 replicates of estimation and prediction results (FLM-SCAD) for all combinations of hyperparameters.

- **cv_summary.R**: This file produces results of Table 1, for rho = 0.75, n = 100, p = 200, q = 10.

