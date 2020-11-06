
if (!require(MRCE)) install.packages('MRCE')
library(MRCE)

i = 1
#args = commandArgs(TRUE)
#eval(parse(text=args))
name0 = 'p120q50_ar1'
name = paste(name0,i,sep='')
n = 100; p = 50; q = 25
Y = read.csv(file=paste("DbHS_",name,"_Y.csv",sep=""),header=F)
X = read.csv(file=paste("DbHS_",name,"_X.csv",sep=""),header=F)
B = read.csv(file=paste("DbHS_",name0,"_B.csv",sep=""),header=F)
Omega_true = read.csv(file=paste("DbHS_",name0,"_Omega.csv",sep=""),header=F)

#the following line needs to be commented when i >= 2
save(B, Omega_true, file = paste("MRCE_",name0,".RData",sep=""))

#for higher dimensions, the range of tuning parameters may need to be made smaller
#for better estimates, use finer intervals
lam1.vec = rev(exp(seq(from=-5, to=0, by=0.5)))
lam2.vec = rev(exp(seq(from=-5, to=0, by=0.5)))
t0 = proc.time()
cvfit = mrce(X=X, Y=Y, lam1.vec = lam1.vec, lam2.vec = lam2.vec, method="cv",
	cov.maxit=30, maxit.out=30, maxit.in=30)
fit = mrce(Y=Y, X=X, lam1=cvfit$best.lam1, lam2=cvfit$best.lam2, method="single",
	cov.maxit=30, maxit.out=30, maxit.in=30)
mrce_time = proc.time()-t0

Bhat = fit$Bhat
Omegahat = fit$omega
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

e = new.env(parent = emptyenv())
load(paste("MRCE_",name0,".RData",sep=""),envir = e)
e$mrce_time[[i]] = mrce_time
e$Bhat[[i]] = Bhat
e$Omegahat[[i]] = Omegahat
e$beta_mse[[i]] = beta_mse
e$nonzerobeta_mse[[i]] = nonzerobeta_mse
e$omega_mse[[i]] = omega_mse
e$diagomega_mse[[i]] = diagomega_mse
e$pred_mse[[i]] = pred_mse
e$avgKL[[i]] = avg_KL

e$SEN_B[[i]] = SEN_B
e$SPE_B[[i]] = SPE_B
e$PREC_B[[i]] = PREC_B
e$ACC_B[[i]] = ACC_B
e$SEN_omega[[i]] = SEN_omega
e$SPE_omega[[i]] = SPE_omega
e$PREC_omega[[i]] = PREC_omega
e$ACC_omega[[i]] = ACC_omega

do.call("save", c(ls(envir = e), list(envir = e, file = paste("MRCE_",name0,".RData",sep=""))))

