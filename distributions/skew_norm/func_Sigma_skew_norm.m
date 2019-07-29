function Sigma = func_Sigma_skew_norm(t,x,delta,xi,n_dim,varargin)
% This function maps the means and variances to \f$\boldsymbol{\Sigma}\f$
% of the multivariate skew normal distribution.
%
% USAGE:
% Sigma = func_Sigma_skewnorm(t,x,xi,n_dim,noise,''additive'')
%
% Parameters:
% t: time vector
% x: vector including means and variances
% xi: (not used, included for consistency and possible extensions)
% n_dim: dimension of measurement
% varargin:
%   * noise: parameter for measurement noise
%   * noisemodel: (so far only ''additive'' supported)
%
% Return values:
% Sigma: (n_t x n_dim x n_dim) \f$\Sigma\f$ of the multivariate skew normal distribution.

noise = zeros(n_dim,1);
noisemodel = 'additive';
if nargin > 5
    noise = varargin{1};
    if nargin > 6
        noisemodel = varargin{2};
    end
end
if ~strcmp(noisemodel,'additive')
    error('Only noise model additive possible!')
end

Sigma = zeros(numel(t),n_dim,n_dim);

D = delta*delta';

for k = 1:numel(t)
    n = n_dim+1;
    for i = 1:n_dim
        for j = i:n_dim
            if isequal(i,j)
                Sigma(k,i,j) = x(k,n)+noise(i)-(1-2/pi)*D(i,j);
            else
                Sigma(k,i,j) = x(k,n)-(1-2/pi)*D(i,j);
                Sigma(k,j,i) = x(k,n)-(1-2/pi)*D(i,j);
            end
            n=n+1;
        end
    end
end

