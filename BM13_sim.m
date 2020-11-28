%-----------------------------------------------------------------------------------------------
% TITLE:        BM13_sim.r reads data, estimates the parameters & calculates some errors by BM13.
%
% VERSION:      1st version (11/20/2020).
%
% AUTHORS:      Yunfan Li, Jyotishka Datta, Bruce A. Craig, Anindya Bhadra,
%
% DESCRIPTION:  This function runs BM13 method for saved simulated data for 50 replicates parallely. For details see:              
%               A. Bhadra, B. K. Mallick, Joint high-dimensional Bayesian variable and covariance selection with an application to eQTL analysis, Biometrics
%               69 (2013) 447–457.
%
% DEPENDS ON:   This depends on a few additional Matlab codes, e.g.
% SSUR_swap.m - these are all saved in the subfolder "helper_files" which
% is added to the search path using the 'addpath' command at the beginning
% of this code. 
%%-----------------------------------------------------------------------------------------------

file_folder = [pwd filesep 'helper_files' filesep];
addpath(file_folder);

i = 1;
name0 = 'p120q50_ar1'
name = [name0 num2str(i)]

data_folder = [pwd filesep 'data' filesep];

n = 100; p = 120; q = 50;

X = csvread(strcat(data_folder,'DbHS_',name,'_X.csv'));
Y = csvread(strcat(data_folder,'DbHS_',name,'_Y.csv'));
B = csvread(strcat(data_folder,'DbHS_',name0,'_B.csv'));
Omega_true = csvread(strcat(data_folder,'DbHS_',name0,'_Omega.csv'));

%%%%%% Set MCMC prior parameters and length of MCMC run%%%%%%%%%%%
%c=0.0005;
mx=X*X';
c=(std(Y(:))/std(mx(:)))^2;
b0=10; D0= .8*eye(q); 
burnin=5000; nmc=10000; 
adj_o=eye(q);

%%% The function below searches for variables and the graph %%%
%%%INPUT: X   = An n X p matrix of predictors
%%%       Y   = An n X q matrix of responses
%%% c,b0,D0   = Prior hyper-parameters
%%%    adj_o  = Initial guess for the graph
%%%burnin,nmc = Burn in and number of MCMC runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%OUTPUT: adj_save   =  Posterior probabilities for the graph (q X q)
%%%        z1_save    =  Posterior probabilities for the variables (p X 1)

t = cputime;
[adj_save,z1_save] = SSUR_swap(Y,X,c,b0,D0,burnin,nmc,adj_o,n,p,q);

%Posterior beta and Sigma will now be regenerated ``conditional''
%on 'z1_o' and 'adj_o' based on a few MCMC runs. This is NOT
%sampling from the full conditional (B | rest) and (Sigma | rest). Instead gamma=gammahat and G = Ghat
%are fixed. So this is empirical Bayes.

%z1post = (z1_save>0.5);
%adjpost = (adj_save>0.5);
%Xgamma = X(:,z1post);
%adj_o = adj_save; z1_o = z1_save(z1post);
Xgamma = X;
adj_o = adj_save; z1_o = z1_save;

%Hyper-parameters
b0=10; D0= .8*eye(q); c = 0.02;

[G]=makedecompgraph(adj_o);

[Lambda, V] =HIWsim(G, b0, D0 , 1); %Initialize Sigma
Bgamma = rMNorm(zeros(q,1), c*V, length(z1_o)); %Initialize Beta
Bgamma = Bgamma';

burnin = 300;
nmc = 500;

for iter =1:nmc+burnin
    
    res = (Y - Xgamma*Bgamma); %find the residual
    S = res'*res;
    
    [Lambda, V] =HIWsim(G, b0+n, D0+S , 1); % posterior of Sigma
    
    postcovar = inv(Xgamma'*Xgamma + (1/c)*eye(length(z1_o)));
    Postcovar = kron(postcovar,inv(Lambda));
    
    Postmean = Postcovar*kron(Xgamma', Lambda)*reshape(Y',q*n,1); 
    mpmean= reshape(Postmean, q,length(z1_o)); 
    mpmean=mpmean';
    
    Bgamma = rMNorm(Postmean, Postcovar,1); %Posterior of beta
    Bgamma = reshape(Bgamma,q, length(z1_o)); %reshape to match dimension
    Bgamma=Bgamma';
   
    if(iter>burnin)
       Bgamma_save(:,:,iter-burnin) = Bgamma;
       V_save(:,:,iter-burnin) = V;
    end
end

BM13_time = cputime-t;

meanB = mean(Bgamma_save,3);
meanV = mean(V_save,3);
Omega_est = inv(meanV);

%MSE of posterior estimates of beta
beta_mse = mean(mean((meanB-B).^2))
nonzerobeta_mse = mean((meanB(B>0)-B(B>0)).^2)
%Squared F_norm of posterior estimates of Omega
omega_mse = mean(mean((Omega_est-Omega_true).^2))
diagomega_mse = mean((diag(Omega_est)-diag(Omega_true)).^2)
%Prediction MSE
test_X = csvread(strcat(data_folder,'DbHS_',name,'_Xtest.csv'));
test_Y = csvread(strcat(data_folder,'DbHS_',name,'_Ytest.csv'));
predict_Y = test_X*meanB;
predict_mse = mean(mean((predict_Y-test_Y).^2))
%Average KL divergence
one = -log(det(omega_mean))+log(det(Omega_true))+trace(omega_mean*inv(Omega_true))-q;
tmp = (X*beta_mean - X*B);
tmp = tmp(:);
two = tmp'*kron(omega_mean,eye(n))*tmp;
avg_KL = one/2 + two/(2*n)

matobj = matfile(strcat('BM13_',name0),'Writable',true);
matobj.BM13_time(1,i) = BM13_time;
matobj.beta_mean(((i-1)*p+1):i*p,1:q) = meanB;
matobj.omega_mean(((i-1)*q+1):i*q,1:q) = Omega_est;
matobj.beta_mse(1,i) = beta_mse;
matobj.nonzerobeta_mse(1,i) = nonzerobeta_mse;
matobj.omega_mse(1,i) = omega_mse;
matobj.diagomega_mse(1,i) = diagomega_mse;
matobj.predict_mse(1,i) = predict_mse;
matobj.avg_KL(1,i) = avg_KL;

%0/1 of posterior \Omega
%omega_zero is a q by q upper triangular matrix, indicating whether posterior prob of inclusion >0.5
omega_zero = true(p,q);
TP_omega = 0; TN_omega = 0; FP_omega = 0; FN_omega = 0;
for l = 1:q
	for j = (l+1):q
		omega_zero(l,j) = adj_save(l,j)<0.5; 
		if Omega_true(l,j)==0 & omega_zero(l,j)==1
			TN_omega = TN_omega+1;
		elseif Omega_true(l,j)==0 & omega_zero(l,j)==0
			FP_omega = FP_omega+1;
		elseif Omega_true(l,j)~=0 & omega_zero(l,j)==1
			FN_omega = FN_omega+1;
		elseif Omega_true(l,j)~=0 & omega_zero(l,j)==0
			TP_omega = TP_omega+1;
		end
	end
end
SEN_omega = TP_omega/(TP_omega+FN_omega)
SPE_omega = TN_omega/(TN_omega+FP_omega)
PREC_omega = TP_omega/(TP_omega+FP_omega)
ACC_omega = (TP_omega+TN_omega)/(TP_omega+TN_omega+FP_omega+FN_omega)		

matobj.SEN_omega(1,i) = SEN_omega;
matobj.SPE_omega(1,i) = SPE_omega;
matobj.PREC_omega(1,i) = PREC_omega;
matobj.ACC_omega(1,i) = ACC_omega;

