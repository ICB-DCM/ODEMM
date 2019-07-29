function drhodxi = func_drhodxi(t,x,dxdxi,xi,varargin)
% This function calculates the derivative of \f$\rho\f$ in case of univariate 
% measurements and a normal distribution assumption.
% 
% USAGE: 
% drhodxi = func_dsigma2dxi_norm(t,x,dxdxi,xi)
%
% Parameters:
% t: time vector (not used, included for consistency and possible extensions)
% x: vector of the means and variances of the observables (not used, 
% included for consistency and possible extensions)
% dxdxi: derivatives of the means and variances of the obsevables
% xi: parameter vector(not used, included for consistency and possible extensions)
%
% Return values:
% drhodxi: derivative of \f$\rho\f$ of the negative binomial distribution.

if nargin > 4
    noise = varargin{1};
    dnoisedxi = varargin{2};
end

if nargin > 4 
    drhodxi = bsxfun(@rdivide,permute(dxdxi(1,:,:),[3,2,1]),noise+x(:,2)) - ...
        bsxfun(@times,bsxfun(@rdivide,x(:,1),(noise+x(:,2)).^2), ... 
        bsxfun(@plus,permute(dxdxi(2,:,:),[3,2,1]),dnoisedxi));    
else
    drhodxi = bsxfun(@rdivide,permute(dxdxi(1,:,:),[3,2,1]),x(:,2)) - ...
        bsxfun(@rdivide,(bsxfun(@times,permute(dxdxi(2,:,:),[3,2,1]),x(:,1))),x(:,2).^2);
end
end