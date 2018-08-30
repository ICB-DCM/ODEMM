function Y = studentstrnd(mu,Sigma,nu,n)
% Random number generation for multivariate t distribution
% according to Pyne (PNAS, 2009).

lambda = gamrnd(nu/2,2/nu,n,1);

for i = 1:n
    Y(:,i) = mvnrnd(mu,Sigma./lambda(i));
end