function f=log_hiwishpdf(X,G,bG,DG)
%LOG_HIWISHPDF Calculate the pdf of HIW(bG,DG) given graph G at the point X

p=size(X,1);
cliques=G{1};
separators=G{2};
numberofcliques=length(cliques);
f=log_iwishpdf(X(cliques(1).ID,cliques(1).ID),bG,DG(cliques(1).ID,cliques(1).ID));

if numberofcliques > 1

for i=2:numberofcliques
    f=f+log_iwishpdf(X(cliques(i).ID,cliques(i).ID),bG,DG(cliques(i).ID,cliques(i).ID))-log_iwishpdf(X(separators(i).ID,separators(i).ID),bG,DG(separators(i).ID,separators(i).ID));
end

elseif numberofcliques==1
    f=f;
end
    
