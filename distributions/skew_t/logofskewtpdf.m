function varargout = logofskewtpdf(varargin)
% Parameters:
%   varargin:
%   y: n_cell x d
%   mu: dx1
%   Sigma: d x d
%   nu: scalar
%   delta: dx1
%   dmudxi: n_xi x d
%   dSigmadxi: n_xi x d x d
%   dnudxi: n_xi x 1
%   ddeltadxi: n_xi x d 

y = varargin{1};
mu = varargin{2};
Sigma = varargin{3};
nu = varargin{4};
delta = varargin{5};
d = size(Sigma,1);
if nargout >= 2
    dmudxi = varargin{6};
    dSigmadxi = varargin{7};
    dnudxi = varargin{8};
    ddeltadxi = varargin{9};
    n_xi = size(dmudxi,1);
end
Omega = Sigma + delta*delta';
invOmega = inv(Omega);
yc = y-mu;
%Z = sum(yc'.*(invOmega*yc)',2);
for i = 1:size(y,2)
    Z(i) = yc(:,i)'*invOmega*yc(:,i);
end
logf = log(2) + logofmvtpdf(y,mu,Omega,nu) + ...
         log(tcdf(((delta'*invOmega*yc)./sqrt(1-delta'*invOmega*delta)).*sqrt((nu+d)./(nu+Z))',nu+d))';
     
varargout{1} = logf;

if nargout>=2
    % n_c x n_xi
    %dZdxi = dZdxi';
    dlogf = nan(n_xi,1);
    for i = 1:n_xi
        dSigmaIndxi = -SigmaIn*squeeze(dSigmadxi(i,:,:))*SigmaIn;
        dZdxi(i) = -yc*SigmaIn*dmudxi(i,:)' - dmudxi(i,:)*SigmaIn*yc' + yc*dSigmaIndxi*yc';
        
        dlogf(i) = 0.5*(bsxfun(@plus,(psi((nu+d)/2)-psi(nu/2)-d/nu)*dnudxi(i) ...
            - trace(SigmaIn*squeeze(dSigmadxi(i,:,:))), ... % n_xi x 1
            -(nu+d)/(nu+Z)*dZdxi(:,i) + (Z*(nu+d)-nu*(nu+Z)*log(1+Z./nu))./(nu*(nu+Z))*dnudxi(i)));
    end
    varargout{2} = dlogf;
end
