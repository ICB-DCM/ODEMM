function dmudxi = func_dmudxi_norm(t,x,dxdxi,Sigma,dSigmadxi,xi,u,n_dim)
% This function calculates the derivative of mu of a normal distribution
%
% Parameters:
% t: time vector (not used, included for consistency and possible extensions)
% x: vector of means and variances (not used, included for consistency and possible extensions)
% dxdxi: derivatives of means and variances
% Sigma: (not used, included for consistency and possible extensions)
% dSigmadxi: (not used, included for consistency and possible extensions)
% xi: parameter vector (not used, included for consistency and possible extensions)
% u: input (not used, included for consistency and possible extensions)
% n_dim: dimension of measurement
%
% Return values:
% dmudxi: derivative of mu of a normal distribution


dmudxi=zeros(numel(t),numel(xi),n_dim);
for k = 1:n_dim
    dmudxi(:,:,k) = permute(dxdxi(k,:,:),[3,2,1]);
end

end
