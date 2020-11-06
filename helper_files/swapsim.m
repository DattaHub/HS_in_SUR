%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Simulated data example              
%%%%%%%%%%%%  for the paper "Joint high-dimensional Bayesian
%%%%%%%%%%%%  variable and covariance selection with an application
%%%%%%%%%%%%  to eQTL analysis", by Bhadra, A. and Mallick, B. K.
%%%%%%%%%%%%  Questions: bhadra@stat.tamu.edu 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
s=RandStream('mt19937ar','seed',123);
RandStream.setDefaultStream(s);

load ./Data/X.txt %X: p=498 X n=120 matrix
X=X';
q=100; % Choose q

trueind=[30 40 57 62 161]; % True predictor variables
truegraph = eye(q);
truegraph(1:30,1:30) =1;   %The true graph, top left 30 X 30 block
Y = simdata(X,q,trueind, truegraph); %Create the simulated data

dimX = size(X);
dimY = size(Y);

n = dimX(1);
p = dimX(2);
q = dimY(2);

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

[adj_save,z1_save] = SSUR_swap(Y,X,c,b0,D0,burnin,nmc,adj_o,n,p,q);

                                                            
savefile = strcat('result.mat');
save(savefile, 'z1_save', 'adj_save');


%%%%%% Plot the result for the variable selection %%%%
set(gca,'fontsize',20);
plot(z1_save, 'linewidth',2);
hold on;
plot(trueind, z1_save(trueind), '.', 'color', 'red', 'MarkerSize',30);
xlabel('Predictor','fontsize',20);
ylabel('Posterior Probability','fontsize',20);
hold off;

%%%%% Plot the result for graph selction %%%%%%
imagesc(adj_save);
colormap(flipud(gray))
colorbar
set(gca,'fontsize',20);
xlabel('Responses','fontsize',20);
ylabel('Responses','fontsize',20);
