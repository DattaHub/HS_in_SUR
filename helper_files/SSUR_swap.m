function [adj_save,z1_save] = SSUR_swap(Y,Xorig,c,b0,D0,burnin,nmc,adj_o,n,porig,q)

rY = reshape(Y',n*q,1);
rY = rY - mean(rY);

priormeanmult=0.4;
phii = 0.01/10;

colindex = [1];
lgamma = length(colindex);
Xnew=zeros(n,lgamma);
rowbegin=1;
colbegin=1;
currcolindex=0;

Xnew = Xorig(:,colindex);
A = eye(n)+c*(Xnew*Xnew');
AL = chol(inv(A));

Z=AL*Y;
Snew = Z'*Z;

nedge = (sum(adj_o(:))-q)/2;
G_o = makedecompgraph(adj_o);

tedge = q*(q-1)/2;
indmx= reshape([1:q^2],q,q); 
upperind=indmx(triu(indmx,1)>0);             
beta=4/(q-1);
w = 1/2;  

adj_save = zeros(q,q);
z1_save = zeros(porig,1);


for iter =1:nmc+burnin
    
    %%%%%%% MH Stochastic search for variables %%%%%
z2=zeros(porig,1);
nvar = length(colindex);
z2(colindex)=1;

if nvar ~= (porig) && nvar > 1 
     delet=(rand(1)<.5);      
     if delet==1  
        varind = find(z2==1);
        delj=randsample(varind,1);
        z2(delj)=0;
        nvar=length(varind);
        prop=nvar/(porig - nvar+1)*(1-phii)/phii;
        else 
        varind = find(z2==0);
        addj=randsample(varind,1);        
        z2(addj)=1;
        nvar=(porig) - length(varind);
        prop=(porig- nvar)/(nvar+1)*(phii)/(1-phii);
        end
        
       elseif nvar==1 | nvar ==0
           addj=randsample(porig,1);
           z2(addj)=1; 
           nvar=nvar+1;
           prop=1/2*(porig-nvar)/(nvar+1)*phii/(1-phii);
       
       elseif nvar==porig   
           delj=randsample(porig,1);
           z2(delj)=0;
           nvar=porig - 1;
           prop=1/2*nvar/(porig - nvar + 1)*(1-phii)/phii;
     end
      
lgamma2=length(find(z2==1));
colindex2=find(z2==1);

Xnew2 = Xorig(:,colindex2);
A2= eye(n) + c*Xnew2*Xnew2';
AL2=chol(inv(A2));
Z2=AL2*Y;
Snew2=Z2'*Z2;

  
logtnewpdf= logtpdf(Snew2,G_o,b0,D0,n);
logtoldpdf=logtpdf(Snew,G_o,b0,D0,n);

rr = logtnewpdf - logtoldpdf + log(prop);
    if rand(1) < exp(rr)
       Snew=Snew2;
       colindex=colindex2;
       lgamma = lgamma2;
       disp 'accept gamma'
       disp 'accept iteration='
       disp(iter)
       disp 'tnewpdf='
       disp (logtnewpdf)
       disp(logtoldpdf)
       disp(log(prop))
       disp(rr)
    end
   

     
 %%%%%%%%%%%%  M-H  stochastic search for graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 adj_star=adj_o;   
                 
 if nedge ~= tedge && nedge ~=0 
     delet=(rand(1)<.5);      
     if delet==1  
        edgeind= upperind(adj_star(upperind)==1);
        delj=edgeind(randsample(nedge,1));
       
        adj_star(delj)=0;
        adj_star=adj_star';
        adj_star(delj)=0;
        prop=nedge/(tedge-nedge+1)*(1-beta)/beta;
        else 
        edgeind= upperind(adj_star(upperind)==0);
        addj=edgeind(randsample((tedge-nedge),1));        
        adj_star(addj)=1;adj_star=adj_star';
        adj_star(addj)=1;
        prop=(tedge-nedge)/(nedge+1)*(beta)/(1-beta);
        end
        
       elseif nedge==0
           edgeind= upperind;
           addj=edgeind(randsample(tedge,1));
           adj_star(addj)=1;adj_star=adj_star';
           adj_star(addj)=1; 
           prop=1/2*(tedge-nedge)/(nedge+1)*beta/(1-beta);
       
       elseif nedge==tedge
           edgeind= upperind;   
           delj=edgeind(randsample(tedge,1));
           adj_star(delj)=0;adj_star=adj_star';
           adj_star(delj)=0;
           prop=1/2*(nedge)/(tedge-nedge+1)*(1-beta)/beta;
     end
      
      if mcard(adj_star)   
      [G_star]=makedecompgraph(adj_star);  
      
      target_star = logtpdf(Snew,G_star,b0,D0,n);
      target_o = logtpdf(Snew, G_o, b0,D0,n);
      
      if rand(1)< (exp(target_star-target_o)*prop)
         adj_o=adj_star;  % Accept proposal
         target_o=target_star;
         G_o=G_star;
         nedge = (sum(adj_o(:))-q)/2;
         disp 'accept G'
         disp 'accept iteration='
         disp(iter)
      end        
      end   
      %%%% end graph%%%%
    
      
      z1=zeros(porig,1);
      z1(colindex)=1;
     %%%%%%%%%%  Below save MCMC moving averages %%%%%%%%%%%%%%
      if(iter>burnin)
          adj_save = ((iter -burnin -1).* adj_save + adj_o)./(iter-burnin);
          z1_save =  ((iter -burnin -1).* z1_save + z1)./(iter-burnin);
      end
end


                                                                                          
