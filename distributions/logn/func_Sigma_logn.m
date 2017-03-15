function Sigma = func_Sigma_logn(t,x,xi,n_dim,varargin)
% This function maps means and variances to \Sigma of a log-normal
% distribution
%
% USAGE:
% Sigma = func_Sigma_logn(t,x,xi,n_dim,noise,noisemodel)
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


% Assign Inputs
noise = zeros(size(x,2),1); % measurement noise
noisemodel = 'multiplicative';
if nargin >= 5
   noise = varargin{1};
   if nargin >=6
       noisemodel = varargin{2};
   end
end
if numel(noise)==1
    disp('Using the same measurement noise for all measurands!');
    noise(2) = noise(1);
end

Sigma = zeros(numel(t),n_dim,n_dim);
switch noisemodel
    case 'additive'
    Sigma(:,1,1) = log((x(:,3)+noise(1))./(x(:,1).^2)+1);
    Sigma(:,2,2) = log((x(:,5)+noise(2))./(x(:,2).^2)+1);
    case 'multiplicative'
    Sigma(:,1,1) = log((x(:,3))./(x(:,1).^2)+1)+noise(1);
    Sigma(:,2,2) = log((x(:,5))./(x(:,2).^2)+1)+noise(2);
    otherwise
    error('Noisemodel not defined!')
end
Sigma(:,1,2) = log(x(:,4)./(x(:,2).*x(:,1))+1);
Sigma(:,2,1) = Sigma(:,1,2);

if ~isreal(Sigma)
    error('Sigma not real!')
end
