function dSigmadxi = func_dSigmadxi_RRE(t,x,dxdxi,xi,n_dim,varargin)
% UNDER CONSTRUCTION!
% This function calculates derivate of Sigma in case of a Reaction Rate
% Equation model
%
% USAGE:
% dSigmadxi = func_dSigmadxi_RRE(t,x,dxdxi,xi,n_dim) \n
% dSigmadxi = func_dSigmadxi_RRE(t,x,dxdxi,xi,n_dim,noise,dnoisedxi)
%
% Parameters: 
% t: time vector
% x: (not used, included for consistency and possible extensions)
% dxdxi: vector including derivatvies of means and variances
% xi: parameter vector
% n_dim: dimension of measurement


dSigmadxi = zeros(numel(t),numel(xi),n_dim,n_dim);
ind = varargin{6};

c=1;
for j = 1:numel(t)
    dSigmadxi(j,:,1,1) = [zeros(1,ind(c)-1),10.^xi(ind(c))*log(10),zeros(1,numel(xi)-ind(c))];
    dSigmadxi(j,:,1,2) = [zeros(1,numel(xi))];
    dSigmadxi(j,:,2,1) = [zeros(1,numel(xi))];
    dSigmadxi(j,:,2,2) = [zeros(1,ind(c+1)-1),10.^xi(ind(c+1))*log(10),zeros(1,numel(xi)-ind(c+1))];
    c=c+2;
end
