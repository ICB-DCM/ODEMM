function dsigma2dxi = func_drhodxi_norm(t,x,dxdxi,xi,varargin)
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

%noisemodel = 'additive';
% if nargin >= 4
%     noise = varargin{1};
%     dnoisedxi = varargin{2};
% end
% if nargin >= 5
%     noisemodel = varargin{3};
% end
% 
% if nargin >= 5 
%     switch noisemodel 
%         case 'additive'
%         dsigma2dxi = bsxfun(@plus,permute(dxdxi(2,:,:),[3,2,1]),dnoisedxi);           
%     end    
% else
    drhodxi = permute(dxdxi(1,:,:),[3,2,1])*rho;
%end
end