function [Zout,varargout] = getLognMeanVar(Z,n_dim,varargin)
% This function calculations the mean and variances, and the corresponding derivatives
% of a multivariate log-normal distribution given the parameters 
% \f$\boldsymbol{\mu}\f$ and \f$\boldsymbol{\Sigma}\f$.
%
% USAGE:
% [Zout] = getLognMeanVar(Z,n_dim) \n
% [Zout,dZdthetaout] = getLognMeanVar(Z,n_dim,dZdtheta)
%
% Parameters:
% Z: (n_dim + n_dim(n_dim+1)/2) x n_t vector with \f$\boldsymbol{\mu}\f$ 
% and  \f$\boldsymbol{\Sigma}\f$
% n_dim: dimension of the multivariate log-normal distribution
% varargin:
% dZdtheta: ((n_dim + n_dim(n_dim+1)/2) x n_theta x n_t derivative)
%
% Return values:
% Zout: (n_dim + n_dim(n_dim+1)/2) x n_t vector with mean and coariance
% dZdthetaout: derivative of Zout

n_t = size(Z,1);
if nargin >=3
    dZdtheta = varargin{1};
    n_theta = size(dZdtheta,2);
end

Zout = nan(size(Z));
mu = Z(:,1:n_dim);
Sigma = nan(n_dim,n_dim,n_t);

n = n_dim+1;
for i = 1:n_dim
    for j = i:n_dim
        Sigma(i,j,:) = Z(:,n);
        n = n+1;
    end
end

if nargout >= 2
    dZdthetaout = nan(size(dZdtheta));
    dmudtheta = dZdtheta(1:n_dim,:,:);
    dSigmadtheta = nan(n_dim,n_dim,n_theta,n_t);
    n = n_dim+1;
    for i = 1:n_dim
        for j = i:n_dim
            dSigmadtheta(i,j,:,:) = dZdtheta(n,:,:);
            n = n+1;
        end
    end
end

for n=1:n_dim
    Zout(:,n) = exp(mu(:,n)+0.5*squeeze(Sigma(n,n,:)));
    if nargout >= 2
        dmutemp = nan(n_theta,n_t);
        dmutemp(:,:) = dmudtheta(n,:,:);
        dSigmatemp = nan(n_theta,n_t);
        dSigmatemp(:,:) = dSigmadtheta(n,n,:,:);
        dZdthetaout(n,:,:) = bsxfun(@times, Zout(:,n)', bsxfun(@plus,dmutemp,0.5*dSigmatemp));
    end
end

n=n_dim+1;
for i = 1:n_dim
    for j = i:n_dim
        
        Sigmai = nan(n_t,1);
        Sigmaj = nan(n_t,1);
        Sigmaij = nan(n_t,1);
        Sigmai(:,1) = Sigma(i,i,:);
        Sigmaj(:,1) = Sigma(j,j,:);
        Sigmaij(:,1) = Sigma(i,j,:);
        Zout(:,n) = exp(mu(:,i)+mu(:,j)+0.5*(Sigmai+Sigmaj)).*(exp(Sigmaij)-1); %n_t x 1
        if nargout >=2
            dmui = nan(n_theta,n_t);
            dmuj = nan(n_theta,n_t);
            dmui(:,:) = dmudtheta(i,:,:);
            dmuj(:,:) = dmudtheta(j,:,:);
            dSigmai = nan(n_theta,n_t);
            dSigmai(:,:) = dSigmadtheta(i,i,:,:);
            dSigmaj = nan(n_theta,n_t);
            dSigmaj(:,:) = dSigmadtheta(j,j,:,:);
            dSigmaij = nan(n_theta,n_t);
            dSigmaij(:,:) = dSigmadtheta(i,j,:,:);
            
            dZdthetaout(n,:,:) = bsxfun(@times,Zout(:,n)',...
                dmui+dmuj+0.5*(dSigmai+dSigmaj))+...
                bsxfun(@times,(exp(mu(:,i)+mu(:,j)+0.5*(Sigmai+Sigmaj)).*exp(Sigmaij))',...
                dSigmaij);
        end
        n=n+1;
    end
end

if nargout>=2
    varargout{1} = dZdthetaout;
end
