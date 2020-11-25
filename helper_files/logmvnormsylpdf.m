function logpdf = logmvnormsylpdf(x,mu,Sigma_var)
% x  = p.n array of n values of p-dim MV normal
% mu = column p vector mean
% Sigma = p.p variance matrix 
% logpdf = n vector of log pdf values
%
[p n]=size(x); 
for i = 1:p 
    e(i)=-((x(i)-mu(i))^2/(2*Sigma_var(i)));
end
logpdf = sum(e) - log(2*pi)*(p/2) - sum(log(sqrt(Sigma_var)));
