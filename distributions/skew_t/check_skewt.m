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
delta = [0.5;0.5];

nu = 3;

%mu = mu(1);
%Sigma = Sigma(1);
%delta = delta(1);
%d=1;
n_xi = 7;

y = mvnrnd(mu,Sigma,n)';


dmudxi = [1,0;
          0,1;
          0,0;
          0,0;
          0,0;
          0,0;
          0,0]';
      
dSigmadxi = zeros(d,d,n_xi);     
dSigmadxi(1,1,3) = 1;        
dSigmadxi(1,2,4) = 1; 
dSigmadxi(2,1,4) = 1;        
dSigmadxi(2,2,5) = 1;   

ddeltadxi = zeros(2,n_xi);
ddeltadxi(1,6) = 1;
ddeltadxi(2,7) = 1;

%logf1 = logofmvnpdf_old(y',mu',Sigma)

%logf2 = logofmvnpdf(y',mu',Sigma)
%% check mvnpdf
% [logf,dlogf] = logofmvnpdf(y,mu,Sigma,dmudxi,dSigmadxi)
% xi = [1,1,1,-0.5,1];
% [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(xi,@(xi) logofmvnpdf(y,[xi(1);xi(2)],[xi(3),xi(4);xi(4),xi(5)],...
%     dmudxi,dSigmadxi,ddeltadxi),1e-4);
% [g,g_fd_f,g_fd_b,g_fd_c]

%%
logf = logofskewtpdf(y,mu,Sigma,nu,delta)
%[logf,dlogf] = logofskewtpdf(y,mu,Sigma,nu,delta,dmudxi,dSigmadxi,dnudxi,ddeltadxi)
%%
% xi = [1,1,1,-0.5,1,1,1];
% [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(xi,@(xi) logofskewnormpdf(y,[xi(1);xi(2)],[xi(3),xi(4);xi(4),xi(5)],[xi(6);xi(7)],...
%     dmudxi,dSigmadxi,ddeltadxi),1e-4);
% [g,g_fd_f,g_fd_b,g_fd_c]

%% generate multivariate skew t random variables
close all
n=1e5;
Omega = Sigma + delta*delta';

w = gamrnd(nu/2,nu/2);
U0 = randn(n,1)./w.^2;
U = mvnrnd(mu,Sigma./w,n)';

for i = 1:n
    Y(:,i) = delta*abs(U0(i)) + U(:,i);
end

%Y = pearsrnd(mu,Sigma,delta,1,n);

%plot(Y(1,:),Y(2,:),'o'); hold on;

meanY = mu + sqrt(2/pi)*delta
mean(Y,2)
covY = Sigma + (1-2/pi)*delta*delta'
cov(Y')
%%
close all
levelsets = [0.001,0.01,0.1];

hs=scatter(Y(1,:),Y(2,:),'.'); hold on;
set(hs,'MarkerEdgeColor','g');
set(hs,'MarkerEdgeAlpha',0.1);
            
[~,kdensity,X1,X2]=kde2d(Y');
contour(X1,X2,kdensity,levelsets,'b'); hold on;

[Y1,Y2] = meshgrid(linspace(min(Y(1,:)),max(Y(1,:)),100),...
    linspace(min(Y(2,:)),max(Y(2,:))),100);
%P = logofmvnpdf([Y1(:)';Y2(:)'],mu,Sigma);
P = exp(logofskewtpdf([Y1(:)';Y2(:)'],mu,Sigma,nu,delta));
P = reshape(P,size(Y1));

contour(Y1,Y2,P,levelsets,'color','r'); hold on;

           

