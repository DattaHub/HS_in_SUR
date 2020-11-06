#The package 'mSSL' has to be installed
#from https://github.com/skdeshpande91/multivariate_SSL
#before running this file.

rm(list = ls())

library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(SSLASSO)
# devtools::install_github(repo = "skdeshpande91/multivariate_SSL/mSSL")

# library(mSSL)
# library(SSLASSO)
library(MASS)
library(mSSL)

base_dir = "C:/Users/jd033/OneDrive/Documents/R/HSGHS"

setwd(paste0(base_dir,"/data"))

# Create cluster for parallel run
library(pacman)
p_load(SSLASSO, MASS, mSSL, parallel, foreach)
no_cores <- detectCores()-1
RNGkind("L'Ecuyer-CMRG")
cl <- makeCluster(no_cores)
# registerDoParallel(cl)
clusterSetRNGStream(cl, 13)

niter = 50

res<- foreach(i = 1:niter, .packages = c("SSLASSO", "MASS","mSSL"))%dopar%{
  name0 = 'p120q50_ar1'
  name = paste0(name0,i)
  n = 100; p = 200; q = 25
  Y = read.csv(file=paste("DbHS_",name,"_Y.csv",sep=""),header=F)
  Y = as.matrix(Y)
  X = read.csv(file=paste("DbHS_",name,"_X.csv",sep=""),header=F)
  X = as.matrix(X)
  B = read.csv(file=paste("DbHS_",name0,"_B.csv",sep=""),header=F)
  B = as.matrix(B)
  Omega_true = read.csv(file=paste0("DbHS_",name0,"_Omega.csv"),header=F)
  Omega_true = as.matrix(Omega_true)

  #the following line needs to be commented when i >= 2
  if(i == 1){
    save(B, Omega_true, file = paste0("CAPME_",name0,".RData"))
  }

  t0 = proc.time()
  #
  #
  fit_mSSL_dcpe <- mSSL_dcpe(X,Y)
  mssl_time = proc.time() - t0

  Bhat = fit_mSSL_dcpe$B
  Omegahat = fit_mSSL_dcpe$Omega
  #MSE of estimated B
  beta_mse = mean((Bhat-as.matrix(B))^2)
  nonzerobeta_mse = mean((Bhat[which(B>0)]-as.matrix(B)[which(B>0)])^2)
  #squared F_norm of estimated Omega_true
  omega_mse = mean((Omegahat-as.matrix(Omega_true))^2)
  diagomega_mse = mean((diag(Omegahat)-diag(as.matrix(Omega_true)))^2)
  #prediction MSE
  testY = read.csv(file=paste("DbHS_",name,"_Ytest.csv",sep=""),header=F)
  testX = read.csv(file=paste("DbHS_",name,"_Xtest.csv",sep=""),header=F)
  pred_Y = as.matrix(testX)%*%Bhat
  pred_mse = mean((pred_Y-as.matrix(testY))^2)
  #average KL divergence
  X = data.matrix(X)
  Y = data.matrix(Y)
  B = data.matrix(B)
  Omega_true = data.matrix(Omega_true)
  one = -log(det(Omegahat))+log(det(Omega_true))+sum(diag((Omegahat%*%solve(Omega_true))))-q;
  tmp = c(X%*%Bhat - X%*%B);
  two = t(tmp)%*%kronecker(Omegahat,diag(n))%*%tmp;
  avg_KL = one/2 + two/(2*n);


  p = dim(X)[2]; q=dim(Y)[2]
  TP_B = 0; TN_B = 0; FP_B = 0; FN_B = 0;
  for (l in 1:p) {
    for (j in 1:q) {
      if (B[l,j]==0 & Bhat[l,j]==0) {
        TN_B = TN_B+1}
      else if (B[l,j]==0 & Bhat[l,j]!=0) {
        FP_B = FP_B+1}
      else if (B[l,j]!=0 & Bhat[l,j]==0) {
        FN_B = FN_B+1}
      else {
        TP_B = TP_B+1}
    }
  }
  SEN_B = TP_B/(TP_B+FN_B)
  SPE_B = TN_B/(TN_B+FP_B)
  PREC_B = TP_B/(TP_B+FP_B)
  ACC_B = (TP_B+TN_B)/(p*q)

  TP_omega = 0; TN_omega = 0; FP_omega = 0; FN_omega = 0;
  for (l in 1:(q-1)) {
    for (j in (l+1):q) {
      if (Omega_true[l,j]==0 & Omegahat[l,j]==0) {
        TN_omega = TN_omega+1}
      else if (Omega_true[l,j]==0 & Omegahat[l,j]!=0) {
        FP_omega = FP_omega+1}
      else if (Omega_true[l,j]!=0 & Omegahat[l,j]==0) {
        FN_omega = FN_omega+1}
      else {
        TP_omega = TP_omega+1}
    }
  }
  SEN_omega = TP_omega/(TP_omega+FN_omega)
  SPE_omega = TN_omega/(TN_omega+FP_omega)
  PREC_omega = TP_omega/(TP_omega+FP_omega)
  ACC_omega = (TP_omega+TN_omega)/(q*(q-1)/2)

  list("mssl_time" = mssl_time, "Bhat" = Bhat, "Omegahat" = Omegahat,
      "beta_mse" = beta_mse, "nonzerobeta_mse" = nonzerobeta_mse, "omega_mse" = omega_mse,
      "diagomega_mse" = diagomega_mse, "pred_mse" = pred_mse, "avg_KL"= avg_KL,
      "SEN_B" = SEN_B, "SPE_B" = SPE_B, "PREC_B" = PREC_B, "ACC_B" = ACC_B,
      "SEN_omega" = SEN_omega, "SPE_omega" = SPE_omega, "PREC_omega"= PREC_omega,
      "ACC_omega" = ACC_omega)
}

stopCluster(cl)

save(res, file = paste0("sim_results_table", name,".RData"))

mssl_results = matrix(0, niter, 11)

for(i in 1:niter){
  mssl_results[i,] = c(res[[i]]$beta_mse, res[[i]]$omega_mse, res[[i]]$pred_mse,
                       +   res[[i]]$avg_KL, res[[i]]$SEN_B, res[[i]]$SPE_B, res[[i]]$PREC_B,
                       +   res[[i]]$SEN_omega, res[[i]]$SPE_omega, res[[i]]$PREC_omega, res[[i]]$mssl_time[3])
}

colnames(mssl_results) = c("beta_mse", "omega_mse", "pred_mse", "avg_KL",
                           "sen_B", "spe_B", "prec_B", "sen_omega", "spe_omega", "prec_omega", "cpu_time")

apply(mssl_results, 2, mean)
apply(mssl_results, 2, sd)

cbind(apply(mssl_results, 2, mean),
apply(mssl_results, 2, sd))
