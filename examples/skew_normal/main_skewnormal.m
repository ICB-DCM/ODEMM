clear all
close all
clc

addpath(genpath('/Users/carolinloos/PhD/ODEMM/ODEMMToolbox/distributions'))
delta = [100;50];
Sigma = [0.5,0.2;0.2,0.3];
eig(Sigma)
mu = [1;2];

Omega = Sigma + delta*delta';
invOmega = inv(Omega);
alpha = delta'*invOmega*(1-delta'*invOmega*delta)^(-1/2);

y = [100;50];
%f = mvnskewnormal(y,mu,Sigma,delta)
logf = logofmvnskewnormal(y,mu,Sigma,delta)
%%
n = 10000;
mu = 1;
Sigma = 1;
delta = 1;
U0 = randn(1,n);
U = zeros(numel(mu),n); 
Omega = Sigma + delta*delta';

for i = 1:n
    U(:,i) = mvnrnd(mu,Omega);
end
Y = bsxfun(@plus,bsxfun(@times,delta,abs(U0)),U);

histogram(Y(1,:),'Normalization','pdf'); hold on;
xgrid = linspace(-100,100,100);
plot(xgrid,exp(logofmvnskewnormal(xgrid,mu,Sigma,delta)));
%%

true_func = @(x) mvnskewnormal(x,mu,Sigma,delta);
x_samples = rand(1,10^6);
sample_value = true_func(x_samples);
max_value = max(sample_value);
accepted = rand(10^6,1) < (sample_value/max_value);
Y = x_samples(accepted);
%%
subplot(3,1,1);
plot(Y(1,:),Y(2,:),'.'); 
subplot(3,1,2);
histogram(Y(1,:));
subplot(3,1,3);
histogram(Y(2,:));
