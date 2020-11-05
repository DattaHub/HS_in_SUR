# HS in SUR

## Matlab and R codes for reproducing the simulation results 

**"Joint mean–covariance estimation via the horseshoe"**
**Authors: Yunfan Li, Jyotishka Datta, Bruce A. Craig, Anindya Bhadra**
Preprint: (https://arxiv.org/abs/1903.06768)[https://arxiv.org/abs/1903.06768]


## Description:

*  DbHS_sim.m: Matlab code for generating data and estimating the parameters $B$ and $\Omega$ by HS_GHS, using the function HSGHS.m. Users can choose between a clique structure or an auto-regressive (AR1) structure for the precision matrix. The `name0` variable should be modified based on the size parameters p and q, and structure 'ar1' or 'clique', e.g. name0 = 'p120q50_clique' or 'p200q25_ar1'. The simulated datasets are automatically saved in a subfolder called '\data\' with an appropriate namestring that can be read by any other program (e.g. BM13_sim.m). This code also generates and saves test data for different error measures such as MSE and prediction error etc. See the main manuscript for descriptions and results. 

* BM13_sim.m reads data, estimates the parameters & calculates some errors by BM13. 

* capme_sim.r reads data, estimates the parameters & calculates some errors by CAPME. Depends on R package CAPME (https://cran.r-project.org/src/contrib/Archive/capme/)

* mrce_sim.r reads data, estimates the parameters & calculates some errors by MRCE. Depends on R package MRCE (https://cran.r-project.org/src/contrib/Archive/MRCE/). 

* mSSL_sim.r reads data, estimates the parameters & calculates some errors by mSSL. Depends on R package mSSL (https://github.com/skdeshpande91/multivariate_SSL). 

* folder 'helper_files' Contains additional files were written by Bhadra and Mallick. They are necessary for running BM13.
