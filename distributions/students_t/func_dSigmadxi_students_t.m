function dSigmadxi = func_dSigmadxi_students_t(t,x,dxdxi,xi,n_dim,varargin)
% This function maps the means and variances to \f$\boldsymbol{\Sigma}\f$
% of the multivariate normal distribution.
%
% USAGE:
% dSigmadxi = func_dSigmadxi_studentst(t,x,dxdxi,xi,n_n_dim)\n
% dSigmadxi = func_dSigmadxi_studentst(t,x,dxdxi,xi,n_n_dim,noise,dnoisedxi,''additive'')

% Parameters:
% t: time vector
% x: (not used, included for consistency and possible extensions)
% dxdxi: vector including derivatvies of the means and variances
% xi: parameter vector
% n_dim: numer of dimensions of multivariate students t distribution
% varargin:
%   * noise: parameter for measurement noise
%   * dnoisedxi: derivative of measurement noise
%   * noisemodel: (so far only ''additive'' supported)
%
% Return values:
% dSigmadxi: (n_t x n_xi x n_dim x n_dim) derivative of
% \f$\boldsymbol{\Sigma}\f$ of the multivariate students t distribution

noise = zeros(n_dim,1);
dnoisedxi = zeros(2,numel(xi));
dSigmadxi = zeros(numel(t),numel(xi),n_dim,n_dim);

noisemodel = 'additive';
if nargin > 5 % measurement noise
    noise = varargin{1};
    dnoisedxi = varargin{2};
    if nargin > 6
        noisemodel = varargin{3};
    end
end
if ~strcmp(noisemodel,'additive')
    error('Only noise_model additive possible!')
end

for k = 1:numel(t)
    n = n_dim+1;
    for i = 1:n_dim
        for j = i:n_dim
            if isequal(i,j)
                dSigmadxi(k,:,i,j) = permute(dxdxi(n,:,k),[3,2,1])+dnoisedxi(i,:);
            else
                dSigmadxi(k,:,i,j) = permute(dxdxi(n,:,k),[3,2,1]);
                dSigmadxi(k,:,j,i) = dSigmadxi(k,:,i,j);
            end
            n=n+1;
        end
    end
end

