function f=log_iwishpdf(X,df,S)


[p,p]=size(X);

    
f=log(det(X))*(-(df+2*p))/2-1/2*trace(inv(X)*S)+log_iwishart_InvA_const(df,S);
