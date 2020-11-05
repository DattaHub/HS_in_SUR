# HS in SUR

## Matlab and R codes for reproducing the simulation results 

### "Joint mean–covariance estimation via the horseshoe"
### Authors: Yunfan Li, Jyotishka Datta, Bruce A. Craig, Anindya Bhadra. 
### Preprint: https://arxiv.org/abs/1903.06768


# Description:

DbHS_sim.m generates data and estimate them by HS_GHS, using the function HSGHS.m.

BM13_sim.m reads data and estimate them by BM13.

capme_sim.r and mrce_sim.r read data and estimate them by CAPME and MRCE.

DbHS_diagnostic.m performs MCMC diagnostic of the sampling.

DbHS_ROC.m, BM13_ROC.m, and DbHS_compare_ROC.m create ROC curves when p=120, q=50, under the clique structured precision matrix.

The other files were written by Bhadra and Mallick. They are necessary for running BM13.
