function Ysim = simdata(X,q,trueind,truegraph)

dimX=size(X);
n = dimX(1);
p = dimX(2);

bsim = 10; Dsim = 0.8*eye(q);
gammasim = trueind;%true variables
csim = 0.3;
adjsim = eye(q);
Xgammasim = X(:,gammasim);


adjsim = truegraph; %true adjacency matrix

if mcard(adjsim)   
    [Gsim]=makedecompgraph(adjsim);  
end

[Lambdasim, Vsim] =HIWsim(Gsim, bsim, Dsim , 1);

m0 = zeros(q,1);
beta = csim.*rMNorm(m0,Vsim, length(gammasim));

Ysim = Xgammasim*beta' + rMNorm(m0,Vsim, n)';
 
 
 