function Y = skewnormrnd(mu,Sigma,delta,n)
% Random number generation for multivariate skew normal distribution
% according to Pyne (PNAS, 2009).

Y = rand(size(Sigma,1),n);
U = mvnrnd(mu,Sigma,n)';
U0 = randn(n,1);

for i = 1:n
    Y(:,i) = delta.*abs(U0(i)) + U(:,i);
end