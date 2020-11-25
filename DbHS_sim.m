%-----------------------------------------------------------------------------------------------
% TITLE:        DbHS_sim.r creates training and test data, estimates the parameters & calculates some errors by HSGHS.
%
% VERSION:      1st version (11/20/2020).
%
% AUTHORS:      Yunfan Li, Jyotishka Datta, Bruce A. Craig, Anindya Bhadra,
%
% DESCRIPTION:  This function runs HSGHS method for saved simulated data for 50 replicates parallely. 
%               For details see the main manuscript.          
%               
% DEPENDS ON:  
%
% DETAILS: 
%       Matlab code for generating data and estimating the parameters $B$ and $\Omega$ by HS_GHS, using the function HSGHS.m.
%       Users can choose between a clique structure or an auto-regressive (AR1) structure for the precision matrix. 
%       The `name0` variable should be modified based on the size parameters p and q, and structure 'ar1' or 'clique', 
%       e.g. name0 = 'p120q50_clique' or 'p200q25_ar1'. The simulated datasets are automatically saved in a subfolder 
%       called '\data\' with an appropriate namestring that can be read by any other program (e.g. BM13_sim.m). 
%       This code also generates and saves test data for different error measures such as MSE and prediction error etc. 
%       See the main manuscript for descriptions and results. 
%------------------------------------------------------------------------------------------------------------------------

%% Generating B and Omega 

clear all

% Uncomment the following lines if you don't want to use a for loop but run
% the code on a server (needs additional set-ups)

%i=getenv('i')
%i=str2num(i)
% i = 1;

data_folder = [pwd filesep 'data' filesep];


for i = 1:50
    name0 = 'p120q50_ar1';
    name = [name0 num2str(i)];

    n = 100; p = 120; q = 50;

    %Random generator for coefficients and precision matrix
    rng(2018)
    %Generate coefficients (p by q)
    %p*q/20=250 elements are nonzero when p=200, q=25
    %p*q/20=300 elements are nonzero when p=120, q=50
    %p*q/5=125 elements are nonzero when p=50, q=25
    %case one: nonzero coefficients~Unif(-2,-0.5) and (0.5,2)
    %case two: all nonzero coefficients equal to five
    
    num_nonzero = p*q/5;
    B = zeros(1,p*q);
    ind = randsample(p*q, num_nonzero);
    B(ind) = 0.5+rand([1 num_nonzero])*1.5;
    B(ind(1:num_nonzero/2)) = -B(ind(1:num_nonzero/2));
    %B(ind) = 5;
    B = reshape(B,p,q);

    %Generate true precision matrix
    
    %True precision matrix
    %%%true omega with an ar1 structure 
    Omega_true = toeplitz([1,0.45,zeros(1,q-2)]);
    %%%true omega with a hubs structure
    %each hub has 10 elements
    %{
    Omega_true = eye(q);
    for a = 1:2
        l = (a-1)*10+1;
        for j = 1:9
            Omega_true(l,l+j) = 0.25;
            Omega_true(l+j,l) = 0.25;
        end
    end
    %}
    %%%true omega with a clique structure; three features in a clique;
    %eight cliques when q=25; nonzero element=0.75
%     {
%     Omega_true = eye(q);
%     ind = randsample(q,24);
%     for a = 1:8
%         clique_ind = ind(((a-1)*3+1):(a*3));
%         for l = 1:3
%             for j = 1:3
%                 Omega_true(clique_ind(l),clique_ind(j)) = 0.75;
%             end
%         end
%     end
%     for a = 1:24
%         Omega_true(ind(a),ind(a)) = 1;
%     end
%     }
    Sigma_true = eye(q)/Omega_true;

    csvwrite(strcat(data_folder,'DbHS_',name0,'_B.csv'),B)  % change to appropriate folder 
    csvwrite(strcat(data_folder,'DbHS_',name0,'_Omega.csv'),Omega_true)
    % Generate X & Y 

    %Random generator for predictors and observed values
    rng(21+i)
    %Generate design matrix (n by p)
    %X has a Toeplitz structure
    tmp_X = zeros(1,p);
    for j = 1:p
        tmp_X(j) = 0.7^(j-1);
    end
    Sigma_X = toeplitz(tmp_X);
    X = mvnrnd(zeros(n,p), Sigma_X);

    %Generate observed values
    E = mvnrnd(zeros(n,q), Sigma_true);
    Y = X*B + E;

    %save data
    csvwrite(strcat(data_folder,'DbHS_',name,'_X.csv'),X) % change to appropriate folder 
    csvwrite(strcat(data_folder,'DbHS_',name,'_Y.csv'),Y)
    
    % create and save test data in the same folder
    rng(48+i)
    m=50;
    test_X = mvnrnd(zeros(m,p), Sigma_X);
    test_E = mvnrnd(zeros(m,q), Sigma_true);
    test_Y = test_X*B + test_E;
    csvwrite(strcat(data_folder,'DbHS_',name,'_Xtest.csv'),test_X) % change to appropriate folder 
    csvwrite(strcat(data_folder,'DbHS_',name,'_Ytest.csv'),test_Y)

end
%% run HS-GHS

%HSHS estimate, burnin=1000 and nmc=5000
burnin = 50; nmc = 200;
t = cputime;
[beta_save,lambda_sq_save,tau_sq_save,...
omega_save,lambda_G_sq_save,tau_G_sq_save] = HSGHS(X,Y,burnin,nmc,eye(q));
DbHS_time = cputime-t;

%% Evaluate on test data 
%%%Create test data
for i = 1:50
    name0 = 'p120q50_ar1'
    name = [name0 num2str(i)];

    rng(48+i)
    m=50;
    test_X = mvnrnd(zeros(m,p), Sigma_X);
    test_E = mvnrnd(zeros(m,q), Sigma_true);
    test_Y = test_X*B + test_E;
    csvwrite(strcat(data_folder,'DbHS_',name,'_Xtest.csv'),test_X) % change to appropriate folder 
    csvwrite(strcat(data_folder,'DbHS_',name,'_Ytest.csv'),test_Y)

end

%%
beta_mean = mean(beta_save,3);
omega_mean = mean(omega_save,3);

%MSE of posterior estimates of beta
beta_mse = mean(mean((beta_mean-B).^2));
nonzerobeta_mse = mean((beta_mean(B>0)-B(B>0)).^2);
%Squared F_norm of posterior estimates of Omega
omega_mse = mean(mean((omega_mean-Omega_true).^2));
diagomega_mse = mean((diag(omega_mean)-diag(Omega_true)).^2);

%Prediction MSE
predict_Y = test_X*beta_mean;
predict_mse = mean(mean((predict_Y-test_Y).^2));

%Average KL divergence
one = -log(det(omega_mean))+log(det(Omega_true))+trace(omega_mean*inv(Omega_true))-q;
tmp = (X*beta_mean - X*B);
tmp = tmp(:);
two = tmp'*kron(omega_mean,eye(n))*tmp;
avg_KL = one/2 + two/(2*n);

matobj = matfile(strcat('DbHS_',name0),'Writable',true);
matobj.DbHS_time(1,i) = DbHS_time;
matobj.beta_mean(((i-1)*p+1):i*p,1:q) = beta_mean;
matobj.omega_mean(((i-1)*q+1):i*q,1:q) = omega_mean;
matobj.beta_mse(1,i) = beta_mse;
matobj.nonzerobeta_mse(1,i) = nonzerobeta_mse;
matobj.omega_mse(1,i) = omega_mse;
matobj.diagomega_mse(1,i) = diagomega_mse;
matobj.predict_mse(1,i) = predict_mse;
matobj.avg_KL(1,i) = avg_KL;

%0/1 of posterior B
%beta_zero is a p by q matrix, indicating whether 50,75,90% CI of beta(i,j) includes zero
k = 1; end_interval = [25 12.5 5];
while k <= 3
	beta_zero = true(p,q);
	TP_B = 0; TN_B = 0; FP_B = 0; FN_B = 0;
	tmp = end_interval(k);
	for l = 1:p
		for j = 1:q
			beta_zero(l,j) = (prctile(beta_save(l,j,:),tmp)<0) & (prctile(beta_save(l,j,:),100-tmp)>0); 
			if B(l,j)==0 & beta_zero(l,j)==1
				TN_B = TN_B+1;
			elseif B(l,j)==0 & beta_zero(l,j)==0
				FP_B = FP_B+1;
			elseif B(l,j)~=0 & beta_zero(l,j)==1
				FN_B = FN_B+1;
			elseif B(l,j)~=0 & beta_zero(l,j)==0
				TP_B = TP_B+1;
			end
		end
	end
	SEN_B(k) = TP_B/(TP_B+FN_B);
	SPE_B(k) = TN_B/(TN_B+FP_B);
	PREC_B(k) = TP_B/(TP_B+FP_B);
	k = k+1;
end
SEN_B
SPE_B
PREC_B

%0/1 of posterior \Omega
%omega_zero is a q by q upper triangular matrix, indicating whether 90% CI of omega(i,j) includes zero
k = 1; end_interval = [25 12.5 5];
while k <= 3
	omega_zero = true(p,q);
	TP_omega = 0; TN_omega = 0; FP_omega = 0; FN_omega = 0;
	tmp = end_interval(k);
	for l = 1:q
		for j = (l+1):q
			omega_zero(l,j) = (prctile(omega_save(l,j,:),tmp)<0) & (prctile(omega_save(l,j,:),100-tmp)>0); 
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
	SEN_omega(k) = TP_omega/(TP_omega+FN_omega);
	SPE_omega(k) = TN_omega/(TN_omega+FP_omega);
	PREC_omega(k) = TP_omega/(TP_omega+FP_omega);
	k = k+1;
end
SEN_omega
SPE_omega
PREC_omega

matobj.SEN_B(i,1:3) = SEN_B;
matobj.SPE_B(i,1:3) = SPE_B;
matobj.PREC_B(i,1:3) = PREC_B;
matobj.SEN_omega(i,1:3) = SEN_omega;
matobj.SPE_omega(i,1:3) = SPE_omega;
matobj.PREC_omega(i,1:3) = PREC_omega;
