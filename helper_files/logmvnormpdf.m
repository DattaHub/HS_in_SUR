function logpdf = logmvnormpdf(x,mu,Sigma)
% x  = p.n array of n values of p-dim MV normal
% mu = column p vector mean
% Sigma = p.p variance matrix 
% logpdf = n vector of log pdf values
%
[p n]=size(x); C=chol(Sigma); e=inv(C)'*(x-repmat(mu,1,n)); 
logpdf = -sum(e.*e)/2 - sum(log(diag(C))) - log(2*pi)*(p/2);
