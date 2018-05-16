clear all
close all
clc
%   y: n_cell x d
%   mu: 1xd
%   Sigma: d x d
%   nu: scalar
%   dmudxi: n_xi x d
%   dSigmadxi: n_xi x d x d
%   dnudxi: n_xi x 1

n = 10;
d = 2;
mu = [1,1]';
Sigma = [1,-0.5;-0.5,1];
nu = 3;


n_xi = 6;
y = mvtrnd(Sigma,nu,n)'+mu;
y = y;
p = logofmvtpdf(y,mu,Sigma,nu)


dmudxi = [1,0;
          0,1;
          0,0;
          0,0;
          0,0;
          0,0]';
      
dSigmadxi = zeros(d,d,n_xi);     
dSigmadxi(1,1,3) = 1;        
dSigmadxi(1,2,4) = 1; 
dSigmadxi(2,1,4) = 1;        
dSigmadxi(2,2,5) = 1;   

dnudxi = zeros(1,n_xi);
dnudxi(6) = 1;

xi = [1,1,1,-0.5,1,3];
[g,g_fd_f,g_fd_b,g_fd_c]=testGradient(xi,@(xi) logofmvtpdf(y,[xi(1);xi(2)],[xi(3),xi(4);xi(4),xi(5)],xi(6),...
    dmudxi,dSigmadxi,dnudxi),1e-4);
[g,g_fd_f,g_fd_b,g_fd_c]


%%
m = mean(y);
C = cov(y)
nu/(nu-2)*Sigma

