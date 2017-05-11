function dsigma2dxi = func_dsigma2dxi_norm(t,x,dxdxi,xi,varargin)
% This function calculates the derivative of \f$\sigma^2\f$ in case of univariate 
% measurements and a normal distribution assumption.
% 
% USAGE: 
% dsigma2dxi = func_dsigma2dxi_norm(t,x,dxdxi,xi)
% dsigma2dxi = func_dsigma2dxi_norm(t,x,dxdxi,xi,noise,dnoisedxi, ''additive'')
%
% Parameters:
% t: time vector (not used, included for consistency and possible extensions)
% x: vector of the means and variances of the observables (not used, 
% included for consistency and possible extensions)
% dxdxi: derivatives of the means and variances of the obsevables
% xi: parameter vector(not used, included for consistency and possible extensions)
% varargin:
%   * noise: parameter for measurement noise
%   * dnoisedxi: derivative of measurement noise
%   * noisemodel: (so far only ''additive'' supported)
%
% Return values:
% dsigma2dxi: derivative of \f$\sigma^2\f$ of the normal distribution.

noisemodel = 'additive';
if nargin > 4
    noise = varargin{1};
    dnoisedxi = varargin{2};
end
if nargin > 5
    noisemodel = varargin{3};
end

if nargin > 5 
    switch noisemodel 
        case 'additive'
        dsigma2dxi = bsxfun(@plus,permute(dxdxi(2,:,:),[3,2,1]),dnoisedxi);           
    end    
else
    dsigma2dxi = permute(dxdxi(2,:,:),[3,2,1]);
end
end