function rt = rMVt(m,df,T)
% rMNorm generates a pxT array  of normal draws
%  from the p-dim N(m,V)
%
p=length(m); 
V=m*m'
rnorm = repmat(reshape(m,p,1),1,T) +chol(V)'*randn(p,T);
rt = rnorm/(chi2rnd(df)/df)


