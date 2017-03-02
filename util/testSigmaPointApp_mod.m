function SP = testSigmaPointApp_mod(varargin)
% Modified version of the testSigmaPointApp_status.m function of the
% SPToolbox.

nonfun = varargin{1};
xi  = varargin{2};
estruct = varargin{3};
op_SP = varargin{4};
op_SP.nderiv = 0;
if(~isfield(op_SP,'nsamples'))
    op_SP.nsamples = 500;
end

[status,SP] = getSigmaPointApp_status_mod(nonfun,xi,estruct,op_SP);

% Dimensions
n_t = size(SP.Y,1);
n_y = size(SP.Y,2);

beta = estruct.beta(xi);
delta = estruct.delta(xi);

[D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,op_SP.type_D);

n_b = size(D,1);

SP.Y_true = NaN(n_t,n_y,op_SP.nsamples);

is = 1;
while is <= op_SP.nsamples
    try
    bsample = mvnrnd(zeros(n_b,1),D);
    %[status,SP.Y_true(:,:,is)] = nonfun(estruct.phi(beta,bsample'));
    [status,~,~,SP.Y_true(:,:,is)] = nonfun(estruct.phi(beta,bsample'));
    is = is + 1
    catch
        is
    end
end

if(any([op_SP.req(1),op_SP.req(1),op_SP.req(4),op_SP.req(5)]))
    SP.my_true = nanmean(SP.Y_true,3);
    DeltaY = bsxfun(@minus,SP.Y_true,SP.my_true);
end

sigma = estruct.sigma_noise(estruct.phi(beta,SP.B_SP(:,1)));

% adapt sigma to proper size
if(op_SP.req(2))
    if(size(sigma,1) == n_t)
        if(size(sigma,2) == 1)
            C_tech = bsxfun(@times,repmat(sigma.^2,[1,n_y,n_y]),permute(eye(n_y),[3,1,2]));
        elseif(size(sigma,2) == n_y)
            C_tech = bsxfun(@times,repmat(sigma.^2,[1,1,n_y]),permute(eye(n_y),[3,1,2]));
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(size(sigma,2) == n_y)
        if(size(sigma,1) == 1)
            C_tech = bsxfun(@times,repmat(sigma.^2,[n_t,1,n_y]),permute(eye(n_y),[3,1,2]));
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(and(size(sigma,1)==1,size(sigma,2)==1))
        C_tech = bsxfun(@times,repmat(sigma.^2,[n_t,n_y,n_y]),permute(eye(n_y),[3,1,2]));
    else
        error('Incompatible size of sigma parametrisation!')
    end
end


% adapt sigma to proper size
if(op_SP.req(5))
    if(size(sigma,1) == n_t)
        if(size(sigma,2) == 1)
            Cz_tech = diag(reshape(repmat(sigma.^2,[1,n_y]),n_t*n_y,1));
        elseif(size(sigma,2) == n_y)
            Cz_tech = diag(reshape(repmat(sigma.^2,[1,1]),n_t*n_y,1));
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(size(sigma,2) == n_y)
        if(size(sigma,1) == 1)
            Cz_tech = diag(reshape(repmat(sigma.^2,[n_t,1]),n_t*n_y,1));
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(and(size(sigma,1)==1,size(sigma,2)==1))
        Cz_tech = diag(reshape(repmat(sigma.^2,[n_t,n_y]),n_t*n_y,1));
    else
        error('Incompatible size of sigma parametrisation!')
    end
end

if(op_SP.req(2))
    SP.Cy_true = 1/op_SP.nsamples*sum(bsxfun(@times,permute(DeltaY,[1,2,4,3]),permute(DeltaY,[1,4,2,3])),4) ...
            + C_tech;
end

if(op_SP.req(3))
    % Cross-covariance
    SP.Cxy_true = 1/op_SP.nsamples*sum(bsxfun(@times,permute(SP.B_SP,[3,1,4,2]),permute(DeltaY,[1,4,2,3])),4);
end

if(op_SP.req(4))
    SP.mz_true = SP.my_true(:);
end

if(op_SP.req(5))
    DeltaZ = reshape(DeltaY,[size(DeltaY,1)*size(DeltaY,2),size(DeltaY,3)]);
    SP.Cz_true = 1/op_SP.nsamples*sum(bsxfun(@times,permute(DeltaZ,[1,3,2]),permute(DeltaZ,[3,1,2])),3) + Cz_tech;  
end

if(op_SP.req(1))
   plotmy(SP.my,SP.my_true,[]) 
end

if(op_SP.req(2))
   plotCy(SP.Cy,SP.Cy_true,[]) 
end

if(op_SP.req(3))
   plotCy(SP.Cxy,SP.Cxy_true,[]) 
end

if(op_SP.req(4))
   plotmz(SP.mz,SP.mz_true,[]) 
end

if(op_SP.req(5))
   plotCz(SP.Cz,SP.Cz_true,[]) 
end



end

