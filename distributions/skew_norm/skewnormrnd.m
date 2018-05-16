function Y = skewnormrnd(mu,Sigma,delta,n)

Omega = Sigma + delta*delta';
D = sqrt(diag(Omega));
invD = inv(D);
Omega_ = invD*Omega*invD;
U = mvnrnd(mu,Omega_,n)';
U0 = randn(n,1);
for i = 1:n
    Y(:,i) = delta.*abs(U0(i)) + U(:,i);
end