function varargout = logofmvtpdf(varargin)
% Modified version of MATLAB function MVNPDF such that the log-density is
% returned. Gradient is added.
%
%MVTPDF multivariate t probability density function (pdf).
%   Y = MVTPDF(X,C,DF) returns the probability density of the multivariate t
%   distribution with correlation parameters C and degrees of freedom DF,
%   evaluated at each row of X.  Rows of the N-by-D matrix X correspond to
%   observations or points, and columns correspond to variables or
%   coordinates.  Y is an N-by-1 vector.
%
%   C is a symmetric, positive definite, D-by-D correlation matrix.  DF is a
%   scalar, or a vector with N elements.
%
%   Note: MVTPDF computes the PDF for the standard multivariate Student's t,
%   centered at the origin, with no scale parameters.  If C is a covariance
%   matrix, i.e. DIAG(C) is not all ones, MVTPDF rescales C to transform it
%   to a correlation matrix.  MVTPDF does not rescale X.
%
%   Example:
%
%      C = [1 .4; .4 1]; df = 2;
%      [X1,X2] = meshgrid(linspace(-2,2,25)', linspace(-2,2,25)');
%      X = [X1(:) X2(:)];
%      p = mvtpdf(X, C, df);
%      surf(X1,X2,reshape(p,25,25));
%
%   See also MVNPDF, MVTCDF, MVTRND, TPDF.

%   Copyright 2005-2011 The MathWorks, Inc.

% Parameters:
%   varargin:
%   y: d x n_cells
%   mu: dx1
%   Sigma: d x d
%   nu: scalar
%   dmudxi: d x n_xi
%   dSigmadxi: d x d x n_xi
%   dnudxi: 1 x n_xi


y = varargin{1};
mu = varargin{2};
Sigma = varargin{3};
nu = varargin{4};

if nargout >= 2
    dmudxi = varargin{5};
    dSigmadxi = varargin{6};
    dnudxi = varargin{7};
    n_xi = size(dmudxi,2);
end

if nargin<4
    error(message('stats:mvtpdf:TooFewInputs'));
elseif ndims(y)~=2
    error(message('stats:mvtpdf:InvalidData'));
end

% Get size of data
y = y';
[d,n] = size(y);
if d<1
    error(message('stats:mvtpdf:TooFewDimensions'));
end

yc=y-mu;
sz = size(Sigma);
if sz(1) ~= sz(2)
    error(message('stats:mvtpdf:BadCorrelationNotSquare'));
elseif ~isequal(sz, [d d])
    error(message('stats:mvtpdf:InputSizeMismatchC'));
end

SigmaIn = inv(Sigma);

logSqrtDetSigma = log(sqrt(det(Sigma)));

Z = sum(yc'.*(SigmaIn*yc)',2);

logNumer = -((nu+d)/2) .* log(1+Z./nu);
logDenom = logSqrtDetSigma + (d/2)*log(nu*pi);
logf = gammaln((nu+d)/2) - gammaln(nu/2) - logDenom+ logNumer;

varargout{1} = logf;

if nargout>=2
    % n_c x n_xi
    %dZdxi = dZdxi';
    dlogf = nan(n_xi,n);
    for i = 1:n_xi
        dSigmaIndxi = -SigmaIn*squeeze(dSigmadxi(:,:,i))*SigmaIn;
        dZdxi(:,i) = -yc'*SigmaIn*dmudxi(:,i) - (dmudxi(:,i)'*SigmaIn*yc)' + sum(yc'.*(dSigmaIndxi*yc)',2);       
        dlogf(i,:) = 0.5*(bsxfun(@plus,(psi((nu+d)/2)-psi(nu/2)-d/nu)*dnudxi(i) ...
            - trace(SigmaIn*squeeze(dSigmadxi(:,:,i))), ... % n_xi x 1
            -(nu+d)./(nu+Z).*dZdxi(:,i) + (Z*(nu+d)-nu*(nu+Z).*log(1+Z./nu))./(nu*(nu+Z))*dnudxi(i)));
    end
    varargout{2} = dlogf;
end
