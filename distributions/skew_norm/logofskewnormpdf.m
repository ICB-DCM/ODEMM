function varargout = logofskewnormpdf(varargin)
% logf: 
% dlogf: 
%
y = varargin{1};
mu = varargin{2};
Sigma = varargin{3};
delta = varargin{4};
d = size(Sigma,1);
if nargout >= 2
    dmudxi = varargin{5};
    dSigmadxi = varargin{6};
    ddeltadxi = varargin{7};
    n_xi = size(dmudxi,2);    
end


n = size(y,1);

Omega = Sigma + delta*delta';
invOmega = inv(Omega); 

a = (1-delta'*invOmega*delta)^(-1/2);
alpha = delta'*invOmega*a;
assert(real(eig(Sigma))>0)

if nargout < 2
    logf = log(2) + logofmvnpdf(y,mu,Omega) + log(normcdf((alpha*(y-mu)')'));
else
    dalphadxi = nan(n_xi,d);
    for i = 1:n_xi
        % to do make more efficient
        dOmegadxi(:,:,i) = dSigmadxi(:,:,i) + ddeltadxi(:,i)*delta'+delta*ddeltadxi(:,i)';
        dinvOmegadxi(:,:,i) = -invOmega*dOmegadxi(:,:,i)*invOmega;
        
        dadxi(i) = 0.5*a*a^(2)*(ddeltadxi(:,i)'*invOmega*delta+...
            delta'*(dinvOmegadxi(:,:,i)*delta + invOmega*ddeltadxi(:,i)));
        
        dalphadxi(i,:) = (ddeltadxi(:,i)'*invOmega + delta'*dinvOmegadxi(:,:,i))*a + ...
            delta'*invOmega*dadxi(i);
    end
    [temp1,dlogofmvnpdfdxi]=logofmvnpdf(y,mu,Omega,dmudxi,dOmegadxi);
    logf = log(2) + temp1 + log(normcdf((alpha*(y-mu)')'));
end


varargout{1}=logf;
if nargout >= 2
    dlogf = nan(n_xi,n);
    for i = 1:n_xi
        dlogf(i,:) = 1./normcdf((alpha*(y-mu)')')'.*normpdf((alpha*(y-mu)')).*...
            (alpha*(-dmudxi(:,i))+(dalphadxi(i,:)*(y-mu))');
    end
    dlogf = dlogf + dlogofmvnpdfdxi; 
    varargout{2}=dlogf;
end