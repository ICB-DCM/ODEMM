function varargout = logofnbinpdf(x,tau,rho,varargin)
% Modified version of MATLAB function MVNPDF such that the log-density is
% returned. Changed parameter names.
%
% USAGE: 
% [q] = logofnbinpdf(x,tau,rho) \\
% [q,dqdxi] = logofnbinpdf(x,tau,rho,dtaudxi,drhodxi)
%

% NBINPDF Negative binomial probability density function.
%   Y = NBINPDF(X,R,P) returns the negative binomial probability density 
%   function with parameters R and P at the values in X.
%   Note that the density function is zero unless X is an integer.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also NBINCDF, NBINFIT, NBININV, NBINLIKE, NBINRND, NBINSTAT, PDF.

%   Copyright 1993-2007 The MathWorks, Inc.


if nargin < 3
    error(message('stats:nbinpdf:TooFewInputs'));
end

if nargin > 3
   dtaudxi = varargin{1};
   drhodxi = varargin{2};
end

[errorcode x tau rho] = distchck(3,x,tau,rho);

if errorcode > 0
    error(message('stats:nbinpdf:InputSizeMismatch'));
end

% Initialize Y to zero.
if isa(x,'single') || isa(tau,'single') || isa(rho,'single')
    q = zeros(size(x),'single');
else
    q = zeros(size(x));
end
if ~isfloat(x)
    x = double(x);
end

% Out of range or missing parameters and missing data return NaN.
% Infinite values for R correspond to a Poisson, but its mean cannot
% be determined from the (R,P) parametrization.
nans = ~(0 < tau & isfinite(tau) & 0 < rho & rho <= 1) | isnan(x);
q(nans) = NaN;

% Negative binomial distribution is defined on the non-negative
% integers.  Data outside this support return 0.
k = find(0 <= x & isfinite(x) & x == round(x)  &  ~nans);
if ~isempty(k)
    lognk = gammaln(tau(k) + x(k)) - gammaln(x(k) + 1) - gammaln(tau(k));
    q(k) = (lognk + tau(k).*log(rho(k)) + x(k).*log1p(-rho(k)));
end
varargout{1}=q;

if nargout == 2
    dqdxi = psi(x(k)+tau(k))*dtaudxi-psi(tau(k))*dtaudxi + log(rho(k))*dtaudxi -...
        x(k)./(1-rho(k))*drhodxi + tau(k)./rho(k)*drhodxi;
    varargout{2}=dqdxi;
end

% 
% % Fix up the degenerate case.
% k = find(p == 1  &  ~nans);
% if ~isempty(k)
%     y(k) = (x(k) == 0);
% end
