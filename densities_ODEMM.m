function varargout = densities_ODEMM(dist,y,varargin)


densities_ODEMM('logn',y,mu{s}(k,:),permute(Sigma{s}(k,:,:),[2,3,1])), ...
    dwdxi{s}(k,n_xi),w{s},k)
switch dist
    case 'norm'
        
    case 'logn'
        mu = varargin{1};
        Sigma  = varargin{2};
        q = bsxfun(@minus,logofmvnpdf(log(y),mu,Sigma),sum(log(y),2));
        if nargout >= 2           
            dmudxi = varargin{3};
            dSigmadxi = varargin{4};
            dwdxi = varargin{5};
            w = varargin{6};
            SigmaIn = inv(Sigma);
            for n_xi = 1:length(xi)
                dSigmaIndxi = -SigmaIn*dSigmadxi(n_xi,:,:)*SigmaIn; 
                H(:,n_xi) = bsxfun(@plus,dwdxi(n_xi), w.*(-0.5)*...
                    (repmat(sum(sum((SigmaIn.').*dSigmadxi(n_xi,:,:))),length(y),1)...
                    +  bsxfun(@minus,mu,log(y))*SigmaIn*dmudxi(n_xi,:)...
                    + (permute(dmudxi(n_xi,:),[1,3,2])*SigmaIn*(bsxfun(@minus,mu,log(y)))')'...
                    + sum((bsxfun(@minus,mu,log(y))*dSigmaIndxi).*bsxfun(@minus,mu,log(y)),2)));
            end
            
        end
    case 't'
        
    case 'neg bin'
        
end


varargout{1} = q;
if nargout >= 2
    varargout{2} = H;
end