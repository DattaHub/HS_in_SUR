# HS in SUR

## Matlab and R codes for reproducing the simulation results 

**"Joint mean–covariance estimation via the horseshoe"**
**Authors: Yunfan Li, Jyotishka Datta, Bruce A. Craig, Anindya Bhadra**
Preprint: [https://arxiv.org/abs/1903.06768](https://arxiv.org/abs/1903.06768) 


## Description:

*  DbHS_sim.m: Matlab code for generating data and estimating the parameters $B$ and $\Omega$ by HS_GHS, using the function HSGHS.m. Users can choose between a clique structure or an auto-regressive (AR1) structure for the precision matrix. The `name0` variable should be modified based on the size parameters p and q, and structure 'ar1' or 'clique', e.g. name0 = 'p120q50_clique' or 'p200q25_ar1'. The simulated datasets are automatically saved in a subfolder called '\data\' with an appropriate namestring that can be read by any other program (e.g. BM13_sim.m). This code also generates and saves test data for different error measures such as MSE and prediction error etc. See the main manuscript for descriptions and results. 

* BM13_sim.m reads data, estimates the parameters & calculates some errors by BM13. For details see:       

          A. Bhadra, B. K. Mallick, Joint high-dimensional Bayesian variable and covariance selection with an application to eQTL analysis, Biometrics, 69 (2013) 447–457.

* capme_sim.r reads data, estimates the parameters & calculates some errors by CAPME. Depends on R package CAPME (https://cran.r-project.org/src/contrib/Archive/capme/). For details see:

         T. T. Cai, H. Li, W. Liu, J. Xie, Covariate-adjusted precision matrix estimation with an application in genetical genomics, Biometrika 100 (2012) 139–156.

* mrce_sim.r reads data, estimates the parameters & calculates some errors by MRCE. Depends on R package MRCE (https://cran.r-project.org/src/contrib/Archive/MRCE/). For details see: 

         A. J. Rothman, E. Levina, J. Zhu, Sparse multivariate regression with covariance estimation, Journal of Computational and Graphical Statistics 19 (2010) 947-962.

* mSSL_sim.r reads data, estimates the parameters & calculates some errors by mSSL. Depends on R package mSSL (https://github.com/skdeshpande91/multivariate_SSL). For details see:
           S. K. Deshpande, V. Rockova, E. I. George, Simultaneous variable and covariance selection with the multivariate spike-and-slab lasso, Journal of Computational and Graphical Statistics 28 (2019) 921-931.

* folder 'helper_files' Contains additional files were written by Bhadra and Mallick. They are necessary for running BM13.
