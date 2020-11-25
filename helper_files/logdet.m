function [a] = logdet(A)

U = chol(A);
a = 2*sum(log(diag(U)));