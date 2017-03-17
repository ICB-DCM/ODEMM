function dmudxi = func_dmudxi_logn_mean(t,x,dxdxi,Sigma,dSigmadxi,xi,u,dim)
% This function calculates the derivative of \f$\boldsymbol{\mu}\f$ of the 
% (multivariate) log-normal distribution for the case of linking the 
% mean of the observables to the mean of the (multivariate) log-normal distribution.
%
% Parameters:
% t: time vector (not used, included for consistency and possible extensions)
% x: vector of the means and covariances (not used, included for consistency and 
% possible extensions)
% dxdxi: derivatives of the means and covariances
% Sigma: (not used, included for consistency and possible extensions)
% dSigmadxi: (not used, included for consistency and possible extensions)
% xi: parameter vector (not used, included for consistency and possible extensions)
% u: input/stimulus (not used, included for consistency and possible extensions)
% n_dim: dimension of measurements
%
% Return values:
% dmudxi: derivative of \f$\boldsymbol{\mu}\f$ of the (multivariate) 
% log-normal distribution

dmudxi=zeros(numel(t),numel(xi),dim);
for k = 1:dim
    dmudxi(:,:,k) = bsxfun(@ldivide,x(:,k),permute(dxdxi(k,:,:),[3,2,1]))-0.5*dSigmadxi(:,:,k,k);
end
end

