function varargout = logofmvnpdf(varargin)
% Modified version of MATLAB function MVNPDF such that the log-density is
% returned.

y = varargin{1};
mu = varargin{2};
Sigma = varargin{3};

if nargout >= 2
    dmudxi = varargin{4};
    dSigmadxi = varargin{5};
    n_xi = size(dmudxi,2);
end


if nargin<1
    error(message('stats:mvnpdf:TooFewInputs'));
elseif ndims(y)~=2
    error(message('stats:mvnpdf:InvalidData'));
end

d = size(Sigma,1);
if ~isequal(size(Sigma,2),d)
    error(message('Sigma has to be a quadratic matrix'))
end

if ~isequal(size(mu,1),d)
    mu = mu';
    y = y';
end

n = size(y,2);
if ~isequal(size(y,1),d)
    error(message('y and mu need to have the same dimension'))
end

% Assume zero mean, data are already centered
if nargin < 2 || isempty(mu)
    yc = y;

% Get scalar mean, and use it to center data
elseif numel(mu) == 1
    yc = y - mu;

% Get vector mean, and use it to center data
elseif ndims(mu) == 2
    [d2,n2] = size(mu);
    if d2 ~= d % has to have same number of coords as X
        error(message('stats:mvnpdf:ColSizeMismatch'));
    elseif n2 == n % lengths match
        yc = y - mu;
    elseif n2 == 1 % mean is a single row, rep it out to match data
        yc = bsxfun(@minus,y,mu);
    elseif n == 1 % data is a single row, rep it out to match mean
        n = n2;
        yc = bsxfun(@minus,y,mu);  
    else % sizes don't match
        error(message('stats:mvnpdf:RowSizeMismatch'));
    end
    
else
    error(message('stats:mvnpdf:BadMu'));
end

% Assume identity covariance, data are already standardized
if nargin < 3 || isempty(Sigma)
    % Special case: if Sigma isn't supplied, then interpret X
    % and Mu as row vectors if they were both column vectors
    if (d == 1) && (numel(X) > 1)
        y = y';
        d = size(y,2);
    end
    xRinv = X0;
    logSqrtDetSigma = 0;
    
% Single covariance matrix
elseif ndims(Sigma) == 2
    sz = size(Sigma);
    if sz(1)==1 && sz(2)>1
        % Just the diagonal of Sigma has been passed in.
        sz(1) = sz(2);
        sigmaIsDiag = true;
    else
        sigmaIsDiag = false;
    end
    
    % Special case: if Sigma is supplied, then use it to try to interpret
    % X and Mu as row vectors if they were both column vectors.
    if (d == 1) && (numel(y) > 1) && (sz(1) == n)
        yc = yc';
        d = size(yc,2);
    end
    
    %Check that sigma is the right size
    if sz(1) ~= sz(2)
        error(message('stats:mvnpdf:BadCovariance'));
    elseif ~isequal(sz, [d d])
        error(message('stats:mvnpdf:CovSizeMismatch'));
    else
        if sigmaIsDiag
            if any(Sigma<=0)
                error(message('stats:mvnpdf:BadDiagSigma'));
            end
            R = sqrt(Sigma);
            xRinv = bsxfun(@rdivide,yc,R);
            logSqrtDetSigma = sum(log(R));
        else
            % Make sure Sigma is a valid covariance matrix
            [R,err] = cholcov(Sigma,0);
            if err ~= 0
                error(message('stats:mvnpdf:BadMatrixSigma'));
            end
            % Create array of standardized data, and compute log(sqrt(det(Sigma)))
            xRinv = yc' / R;
            logSqrtDetSigma = sum(log(diag(R)));
        end
    end
end   


% The quadratic form is the inner products of the standardized data
quadform = sum(xRinv.^2, 2);

logf = -0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2;

varargout{1}=logf;

if nargout >= 2
    dlogf = nan(n_xi,n);
    SigmaIn = inv(Sigma);
    for i = 1:n_xi
        dSigmaIndxi = -SigmaIn*squeeze(dSigmadxi(:,:,i))*SigmaIn;
        dZdxi(:,i) = -yc'*SigmaIn*dmudxi(:,i) - (dmudxi(:,i)'*SigmaIn*yc)' + sum(yc'.*(dSigmaIndxi*yc)',2);       
        
        dlogf(i,:) = -0.5*(trace(SigmaIn*squeeze(dSigmadxi(:,:,i)))+dZdxi(:,i))';
    end
    varargout{2} = dlogf;
end