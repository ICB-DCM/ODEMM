function dmudxi = func_dmudxi_skew_norm(t,x,dxdxi,delta,ddeltadxi,xi,u,n_dim)
% This function calculates the derivative of \f$\boldsymbol{\mu}\f$ of the
% (multivariate) skew normal distribution.
%
% Parameters:
% t: time vector (not used, included for consistency and possible extensions)
% x: vector of means and variances (not used, included for consistency and possible extensions)
% dxdxi: derivatives of means and variances
% delta: 
% ddeltadxi:
% xi: parameter vector (not used, included for consistency and possible extensions)
% u: input (not used, included for consistency and possible extensions)
% n_dim: dimension of measurement
%
% Return values:
% dmudxi: derivative of \f$\boldsymbol{\mu}\f$ of the (multivariate) 
% skew normal distribution.

dmudxi=zeros(numel(t),numel(xi),n_dim);
for k = 1:n_dim
    dmudxi(:,:,k) = permute(dxdxi(k,:,:),[3,2,1]) - sqrt(2/pi)*ddeltadxi(k,:);
end

end
