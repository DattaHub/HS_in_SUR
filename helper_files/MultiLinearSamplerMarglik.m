function [theta_save,C_save,Lambda_save,V_save,logm1,logm2] = MultiLinearSamplerMarglik(Y,Fnp,Fpn,n,p,T,m0,nu0,nu1,b0,D0,burnin,nmc,z1,adj)

G = makedecompgraph(adj);
C0_var = nu0.^2;
C0_var(z1) = nu1(z1).^2;
C0 = diag(C0_var);




% n = length(z1);
F1 = zeros(n,T*p);


% % %%%%  MLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F1 = Fnp;
% invC = F1*Fnp';
% C = inv(invC);
% m_MLE =  C* (F1* reshape(Y,p*T,1) );
% % m([1,6,11,16])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = zeros(n,1); % Initial value 
% 
% theta(1:9) = result(1).beta; % Initial value 
% theta(10:18) = result(2).beta; % Initial value 
% theta(19:27) = result(3).beta; % Initial value 
% theta(28:36) = result(4).beta; % Initial value 
% 
% theta = theta([9,1:2:8,2:2:8,18,10:2:17,11:2:17,27,19:2:26,20:2:26, 36,28:2:35,29:2:35]);
% % 

theta_save = zeros(n,nmc);
Lambda_save = zeros(p,p,nmc);
V_save = Lambda_save;
m_save  = zeros(n,nmc) ;
C_save  = zeros(n,n,nmc);
DT_save = zeros(p,p,nmc);
       
       
       
for iter =1:nmc+burnin

e = Y - reshape(Fnp'*theta,p,T);
DT = D0 + e*e';
[Lambda,V] = HIWsim(G,b0+T,DT,1);

% Lambda = wishart_InvA_rnd(b0+T+p-1,DT,1);

F10 = Fpn'*Lambda;
for i=1:T
F1(:,(i-1)*p+1:i*p) = F10((i-1)*n+1:i*n,:);
end

invC = F1*Fnp' + inv(C0);
% invC = sum_prod(F,Lambda)+inv(C0);

C = inv(invC);
C = (C+C')/2;
m =  C* (inv(C0)*m0 + F1* reshape(Y,p*T,1) );
theta = rMNorm(m,C,1);
     
   if(iter>burnin)
       theta_save(:,iter-burnin) = theta;
       Lambda_save(:,:,iter-burnin) = Lambda;
       V_save(:,:,iter-burnin) = V;
       m_save(:,iter-burnin)= m ;
       C_save(:,:,iter-burnin)= C;
       DT_save(:,:,iter-burnin) = DT;
   end

end % end for iter =1:nmc+burnin





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Chibs methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lambda_plugin = mean(Lambda_save,3);
V_plugin = inv(Lambda_plugin);

% F1 = Fnp * kron(eye(T),Lambda_plugin');

F10 = Fpn'*Lambda_plugin;
for i=1:T
F1(:,(i-1)*p+1:i*p) = F10((i-1)*n+1:i*n,:);
end

invC = F1*Fnp' + inv(C0);



invC = F1*Fnp' + inv(C0);

C_plugin = inv(invC);
m_plugin  =  C_plugin* (inv(C0)*m0 + F1* reshape(Y,p*(T),1) );

vecY = reshape(Y,p*(T),1);


lognumerator1 =  -(p*(T))/2*log(2*pi)+(T)/2*logdet(Lambda_plugin)-1/2*logdet(C0)+1/2*logdet(C_plugin) ...
    -1/2* reshape(Lambda_plugin*Y,1,T*p)*vecY  -1/2*m0'*inv(C0)*m0  +1/2*m_plugin'*inv(C_plugin)*m_plugin ...
    + log_hiwishpdf(V_plugin,G,b0,D0);


theta_plugin = mean(theta_save,2);
e = Y - reshape(Fnp'*theta_plugin,p,T);
lognumerator2 =  -(p*(T))/2*log(2*pi)+log_hiwishart_InvA_const(G,b0,D0)-log_hiwishart_InvA_const(G,b0+T,D0+e*e') ...
    + logmvnormpdf(theta_plugin,m0,C0);


logLambda = zeros(1,nmc);
logTheta = zeros(1,nmc);

for iter =1:nmc
 logLambda(iter) = log_hiwishpdf(V_plugin,G,b0+T,DT_save(:,:,iter));
 logTheta(iter) = logmvnormpdf(theta_plugin,m_save(:,iter),C_save(:,:,iter));
end

logdenominator1  =  log(sum(exp(logLambda-max(logLambda)))) +  max(logLambda)-log(nmc);
logdenominator2  =  log(sum(exp(logTheta-max(logTheta)))) +  max(logTheta)-log(nmc);

logm1 = lognumerator1 - logdenominator1;
logm2 = lognumerator2 - logdenominator2;








