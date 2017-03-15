function Sigma = func_Sigma_norm(t,x,xi,n_dim,varargin)
% This function maps the means and variances to Sigma of a normal
% distribution
%
% USAGE:
% Sigma = func_Sigma_norm(t,x,xi,n_dim,noise,noisemodel)
%
% Parameters:
% t: time vector
% x: vector including means and variances
% xi: (not used, included for consistency and possible extensions)
% n_dim: dimension of measurement
% varargin: 
%   noise:
%   noisemodel:
%
% Return values:
% Sigma: n_t x n_dim x n_dim 

noise = [0;0];
noisemodel = 'additive';
if nargin >= 4
   noise = varargin{1};
   if nargin >=5
       noisemodel = varargin{2};
   end
end
if ~strcmp(noisemodel,'additive')
    error('Only noisemodel additive possible!')
end
Sigma = zeros(numel(t),n_dim,n_dim);
if ~isequal(n_dim,2)
    error('Todo: adapt for n_dim > 2')
end
if n_dim==2
    for k=1:numel(t)
     Sigma(k,:,:) = [x(k,3)+noise(1),x(k,4); x(k,4),x(k,5)+noise(2)];
    end
else
    for k = 1:numel(t)
        n = n_dim+1;
        for i = 1:n_dim
            for j = 1:n_dim
                if isequal(i,j)
                    Sigma(k,i,j) = x(k,n)+noise(i);
                else                   
                    Sigma(k,i,j) = x(k,n);
                    Sigma(k,j,i) = x(k,n);
                end
                n=n+1;
            end
        end
    end
end
