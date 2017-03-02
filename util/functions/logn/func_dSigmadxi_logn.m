function dSigmadxi = func_dSigmadxi_logn(t,x,dxdxi,xi,n_dim,varargin)
% This function maps means and variances to \Sigma of a log-normal
% distribution
%
% USAGE:
% dSigmadxi = func_dSigmadxi_logn(t,x,dxdxi,xi,n_n_dim)
% dSigmadxi = func_dSigmadxi_logn(t,x,dxdxi,xi,n_n_dim,noise,dnoisedxi,'multiplicative')

% Parameters:
% t: time vector
% x: vector including means and variances
% dxdxi: vector including derivatvies of means and variances
% xi: parameter vector
% n_dim: numer of dimensions of multivariate log-normal distribution
% varargin:
%   noise:
%   dnoisedxi:
%   noisemodel:
%       = ''multiplicative''
%       = ''additive''
%
% Return values:
% dSigmadxi: n_t x n_xi x n_dim x n_dim


noise = zeros(size(x,2),1);
dnoisedxi = zeros(size(x,2),numel(xi));
noisemodel = 'multiplicative';

if nargin >= 6 % measurement noise
    noise = varargin{1};
    dnoisedxi = varargin{2};
    if nargin >=8
        noisemodel = varargin{3};
    end
end

if numel(noise)==1
    disp('Using the same measurement noise for all measurands!');
    noise(2) = noise(1);
    dnoisedxi = [dnoisedxi,dnoisedxi]';
end

if ~isequal(n_dim,2)
    error('TODO: adapt for n_dim > 2');
end
dSigmadxi = zeros(numel(t),numel(xi),n_dim,n_dim);

for k = 1:numel(t)
    m1 = x(k,1);
    m2 = x(k,2);
    C11 = x(k,3);
    C12 = x(k,4);
    C22 = x(k,5);
    switch noisemodel
        case 'additive'
            for j = 1:numel(xi)
                dSigmadxi(k,j,1,1) = 1/(((C11+noise(1))/m1^2)+1)*(m1^2*(permute(dxdxi(3,j,k),[3,2,1])+dnoisedxi(1,j))...
                    -(C11+noise(1))*(m1*permute(dxdxi(1,j,k),[3,2,1])+m1*permute(dxdxi(1,j,k),[3,2,1])))/(m1^4);
                dSigmadxi(k,j,2,2) = 1/(((C22+noise(2))/m2^2)+1)*(m2^2*(permute(dxdxi(5,j,k),[3,2,1])+dnoisedxi(2,j))...
                    -(C22+noise(2))*(m2*permute(dxdxi(2,j,k),[3,2,1])+m2*permute(dxdxi(2,j,k),[3,2,1])))/(m2^4);
                dSigmadxi(k,j,1,2) = 1/((C12/(m1*m2))+1)*(m1*m2*permute(dxdxi(4,j,k),[3,2,1])...
                    -C12*(m2*permute(dxdxi(1,j,k),[3,2,1])+m1*permute(dxdxi(2,j,k),[3,2,1])))/((m1*m2)^2);
                dSigmadxi(k,j,2,1) = dSigmadxi(k,j,1,2);
            end
        case 'multiplicative'
            for j = 1:numel(xi)
                dSigmadxi(k,j,1,1) = 1/((C11/m1^2)+1)*(m1^2*permute(dxdxi(3,j,k),[3,2,1])...
                    -C11*(m1*permute(dxdxi(1,j,k),[3,2,1])+m1*permute(dxdxi(1,j,k),[3,2,1])))/(m1^4) + dnoisedxi(1,j);
                dSigmadxi(k,j,2,2) = 1/((C22/m2^2)+1)*(m2^2*permute(dxdxi(5,j,k),[3,2,1])...
                    -C22*(m2*permute(dxdxi(2,j,k),[3,2,1])+m2*permute(dxdxi(2,j,k),[3,2,1])))/(m2^4) +dnoisedxi(2,j);
                dSigmadxi(k,j,1,2) = 1/((C12/(m1*m2))+1)*(m1*m2*permute(dxdxi(4,j,k),[3,2,1])...
                    -C12*(m2*permute(dxdxi(1,j,k),[3,2,1])+m1*permute(dxdxi(2,j,k),[3,2,1])))/((m1*m2)^2);
                dSigmadxi(k,j,2,1) = dSigmadxi(k,j,1,2);
            end
        otherwise
            error('Noisemodel not defined!')
    end
    
end
end

