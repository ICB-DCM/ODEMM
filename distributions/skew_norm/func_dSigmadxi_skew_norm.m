function dSigmadxi = func_dSigmadxi_skew_norm(t,x,dxdxi,delta,ddeltadxi,xi,n_dim,varargin)
% This function maps the means and variances to \f$\boldsymbol{\Sigma}\f$
% of the multivariate skew normal distribution.
%
% USAGE:
% dSigmadxi = func_dSigmadxi_skewnorm(t,x,dxdxi,xi,n_n_dim)\n
% dSigmadxi = func_dSigmadxi_skewnorm(t,x,dxdxi,xi,n_n_dim,noise,dnoisedxi,''additive'')

% Parameters:
% t: time vector
% x: (not used, included for consistency and possible extensions)
% dxdxi: vector including derivatives of the means and variances
% xi: parameter vector
% n_dim: numer of dimensions of multivariate skew normal distribution
% varargin:
%   * noise: parameter for measurement noise
%   * dnoisedxi: derivative of measurement noise
%   * noisemodel: (so far only ''additive'' supported)
%
% Return values:
% dSigmadxi: (n_t x n_xi x n_dim x n_dim) derivative of
% \f$\boldsymbol{\Sigma}\f$ of the multivariate skew normal distribution

noise = zeros(n_dim,1);
dnoisedxi = zeros(n_dim,numel(xi));
dSigmadxi = zeros(numel(t),numel(xi),n_dim,n_dim);

noisemodel = 'additive';
if nargin >= 8 % measurement noise
    noise = varargin{1};
    dnoisedxi = varargin{2};
    if nargin >= 9
        noisemodel = varargin{3};
    end
end
if ~strcmp(noisemodel,'additive')
    error('Only noisemodel additive possible!')
end


dDdxi = zeros(numel(xi),n_dim,n_dim);
for i = 1:numel(xi)
    dDdxi(i,:,:) = ddeltadxi(:,i)*delta'+ delta*ddeltadxi(:,i)';
end

for k = 1:numel(t)
    n = n_dim+1;
    for i = 1:n_dim
        for j = i:n_dim
            if isequal(i,j)
                dSigmadxi(k,:,i,j) = permute(dxdxi(n,:,k),[3,2,1])+dnoisedxi(i,:)-(1-2/pi)*dDdxi(:,i,j)';
            else
                dSigmadxi(k,:,i,j) = permute(dxdxi(n,:,k),[3,2,1])-(1-2/pi)*dDdxi(:,i,j)';
                dSigmadxi(k,:,j,i) = permute(dxdxi(n,:,k),[3,2,1])-(1-2/pi)*dDdxi(:,j,i)';
            end
            n=n+1;
        end
    end
end

