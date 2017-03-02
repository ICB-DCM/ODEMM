function dsigma2dxi = func_dsigma2dxi_norm(t,x,dxdxi,xi,varargin)
% This function calculates the derivative of sigma^2 in case of univariate 
% measurements and a normal distribution assumption
% 
% USAGE: 
% dsigma2dxi = func_dsigma2dxi_norm(t,x,dxdxi,xi)
% dsigma2dxi = func_dsigma2dxi_norm(t,x,dxdxi,xi,noise,dnoisedxi,noise_model)
%
% Parameters:
% t: time vector (not used, included for consistency and possible extensions)
% x: vector of means and variances (not used, included for consistency and possible extensions)
% dxdxi: derivatives of means and variances
% xi: parameter vector(not used, included for consistency and possible extensions)
% varargin:
%   noise: parameter for measurement noise
%   dnoisedxi: derivative of measurement noise
%   noise_model: 
%
% Return values:
% dsigma2dxi: derivative of sigma^2 of a normal distribution

noise_model = 'additive';
if nargin >= 4
    noise = varargin{1};
    dnoisedxi = varargin{2};
end
if nargin >= 5
    noise_model = varargin{3};
end

if nargin >= 5 
    switch noise_model 
        case 'additive'
        dsigma2dxi = bsxfun(@plus,permute(dxdxi(2,:,:),[3,2,1]),dnoisedxi);           
    end    
else
    dsigma2dxi = permute(dxdxi(2,:,:),[3,2,1]);
end
end