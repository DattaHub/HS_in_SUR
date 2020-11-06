function f=log_tpdf(X,df,S,n)


[p,p]=size(X);


f=1/2*(df+p-1)*log(det(S)) -1/2*(df+p+n-1)*log(det(X+S)) + ...
  log_multi_gamma(p,(df+p+n-1)/2) - log_multi_gamma(p,(df+p-1)/2) ...
  -(p*n/2)*log(pi);





