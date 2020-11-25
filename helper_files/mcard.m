function [isdecomp] = mcard(A)
% maximum cardinality search of a graph given adjmx A on p nodes; here
%   we make sure to set the diag of A to 1's as per our convention
p=size(A,1); A=A-diag(diag(A))+eye(p);
isdecomp=1;
Nnodes=1; nnum=1; UNnodes=2:p;           % choose to number node 1 to start, arbitrarily
while (isdecomp & nnum<p)
   % find numbers of numbered neighbours of unnumbered nodes
   nn = sum(A(UNnodes,Nnodes),2);
   [a,j]=sort(nn,'descend');  j=j(1);
   Nnodes=[Nnodes UNnodes(j)]; UNnodes(j)=[]; j=Nnodes(end);
   nej=intersect(Nnodes,find(A(j,:)));
   isdecomp = all(all(A(nej,nej)));
   nnum=nnum+1;
end
