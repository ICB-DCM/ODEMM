function Sigma = func_Sigma_RRE(varargin)
% UNDER CONSTRUCTION!
t = varargin{1};
x = varargin{2};
xi = varargin{3};
n_dim = varargin{4};
ind = varargin{5};
noise = 0;
if nargin >= 6
   noise = varargin{6};
end
Sigma = zeros(numel(t),n_dim,n_dim);
c=1;
for j=1:numel(t)
  Sigma(j,:,:) = [10.^xi(ind(c)),0;0,10.^xi(ind(c+1))];
  c=c+2;
end
