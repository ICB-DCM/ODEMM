function dsigma2dxi = func_dsigma2dxi_logn(t,x,dxdxi,xi,varargin)
% This function calcuatates the derivative of sigma^2 in case of univariate
% measurements and a log-normal distribution assumption
% 
% USAGE: 
% dsigma2dxi = func_dsigma2dxi_logn(t,x,dxdxi,xi)
% dsigma2dxi = func_dsigma2dxi_logn(t,x,dxdxi,xi,noise,dnoisedxi,noise_model)
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
% dsigma2dxi: derivative of sigma^2 of a log-normal distribution

noise_model = 'multiplicative';
if nargin >= 5
    noise = varargin{1};
    dnoisedxi = varargin{2};
end
if nargin >= 6
    noise_model = varargin{3};
end

if nargin >= 5
    switch noise_model 
        case 'multiplicative'
        dsigma2dxi = bsxfun(@plus,bsxfun(@times, 1./(x(:,2)./(x(:,1).^2)+1),...
        bsxfun(@rdivide, bsxfun(@times,(x(:,1).^2), permute(dxdxi(2,:,:),[3,2,1])) -...
        bsxfun(@times, x(:,2).*2.*x(:,1), permute(dxdxi(1,:,:),[3,2,1])), x(:,1).^4 )),dnoisedxi);
        case 'additive'
        dsigma2dxi = bsxfun(@times, 1./((x(:,2)+noise)./(x(:,1).^2)+1),...
        bsxfun(@rdivide, bsxfun(@times,(x(:,1).^2), bsxfun(@plus,permute(dxdxi(2,:,:),[3,2,1]),dnoisedxi)) -...
        bsxfun(@times, (x(:,2)+noise).*2.*x(:,1), permute(dxdxi(1,:,:),[3,2,1])), x(:,1).^4 ));           
    end    
else
    dsigma2dxi = bsxfun(@times, 1./(x(:,2)./(x(:,1).^2)+1),...
    bsxfun(@rdivide, bsxfun(@times,(x(:,1).^2), permute(dxdxi(2,:,:),[3,2,1])) -...
    bsxfun(@times, x(:,2).*2.*x(:,1), permute(dxdxi(1,:,:),[3,2,1])), x(:,1).^4 ));
end
end