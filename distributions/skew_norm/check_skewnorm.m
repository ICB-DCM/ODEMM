clear all
close all
clc

d = 2;
mu = [1,1];

Sigma = [3,-0.5;-0.5,1];
delta = [-20;20];

d=2;
n_xi = 7;
n=1;

y = mvnrnd(mu,Sigma,n);

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


%%
logf = logofskewnormpdf(y,mu,Sigma,delta)
[logf,dlogf] = logofskewnormpdf(y,mu,Sigma,delta,dmudxi,dSigmadxi,ddeltadxi)
%%
%xi = [1,1,1,-0.5,1,1,1];
xi = rand(8,1);
xi(3:5) = [1,-0.5,1];
[g,g_fd_f,g_fd_b,g_fd_c]=testGradient(xi,@(xi) logofskewnormpdf(y,[xi(1),xi(2)],[xi(3),xi(4);xi(4),xi(5)],[xi(6);xi(7)],...
    dmudxi,dSigmadxi,ddeltadxi),1e-4);
[g,g_fd_f,g_fd_b,g_fd_c]

%% check 1D
d = 1;
n_xi = 3;
n=1;
mu = 200;
Sigma = 300;
delta = -100;
y = mvnrnd(mu,Sigma,n);
dmudxi = [1;
          0;
          0]';
            
y = 444;
mu = 490;
Sigma = 1.2e4;
delta = -13;
      
dSigmadxi = zeros(d,d,n_xi);     
dSigmadxi(1,1,2) = 1;  

ddeltadxi = zeros(1,n_xi);
ddeltadxi(1,3) = 1;


%%
logf = logofskewnormpdf(y,mu,Sigma,delta)
[logf,dlogf] = logofskewnormpdf(y,mu,Sigma,delta,dmudxi,dSigmadxi,ddeltadxi)

xi = rand(3,1);
[g,g_fd_f,g_fd_b,g_fd_c]=testGradient(xi,@(xi) logofskewnormpdf(y,[xi(1)],[xi(2)],[xi(3)],...
    dmudxi,dSigmadxi,ddeltadxi),1e-4);
[g,g_fd_f,g_fd_b,g_fd_c]



%% generate multivariate skew normal random variables
close all
n=1e6;
mu = [10;0];
delta = [20;10];
Sigma = eye(2)*2;

Y = skewnormrnd(mu,Sigma,delta,n);
meanY = mu + sqrt(2/pi)*delta
mean(Y,2)

covY = Sigma + (1-2/pi)*delta*delta'
cov(Y')
%%
subplot(1,3,1);
scatter(Y(1,:),Y(2,:)); axis square;

histogram(Y,'normalization','pdf'); hold on;
grid = linspace(min(Y),max(Y),1000);
p = exp(logofskewnormpdf(grid,mu,Sigma,delta));
plot(grid,exp(logofskewnormpdf(grid,mu,Sigma,delta)),'LineWidth',2);

%xlim([0,100])
%%
f = @(y) exp(logofskewnormpdf(y,mu,Sigma,delta))
sampleDist(f,1e2,1e2,[-1e2,1e2]) 
%%
close all
levelsets = 5;%[0.001,0.01,0.1];

hs=scatter(Y(1,:),Y(2,:),'.'); hold on;
set(hs,'MarkerEdgeColor','g');
set(hs,'MarkerEdgeAlpha',0.1);
            
[~,kdensity,X1,X2]=kde2d(Y');
contour(X1,X2,kdensity,levelsets,'b'); hold on;

[Y1,Y2] = meshgrid(linspace(min(Y(1,:)),max(Y(1,:)),100),...
    linspace(min(Y(2,:)),max(Y(2,:))),100);
%P = logofmvnpdf([Y1(:)';Y2(:)'],mu,Sigma);
%%
P = exp(logofskewnormpdf([Y1(:)';Y2(:)'],mu,Sigma,delta));
P = reshape(P,size(Y1));

contour(Y1,Y2,P,levelsets,'color','r'); hold on;
%%

           
%%
%% generate multivariate skew normal random variables
close all
n=1e6;

Y = mvnrnd(mu,Sigma,n)';

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
P = exp(logofmvnpdf([Y1(:)';Y2(:)'],mu,Sigma,delta));
P = reshape(P,size(Y1));

contour(Y1,Y2,P,levelsets,'color','r'); hold on;

           

%%
close all
grid = linspace(0,150,100);
plot(grid,exp(logofskewnormpdf(grid',10,10,10)),'LineWidth',2); hold on;
plot(grid,exp(logofskewnormpdf(grid',40,10,20)),'LineWidth',2);
plot(grid,exp(logofskewnormpdf(grid',100,10,30)),'LineWidth',2);

