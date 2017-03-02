function [status,SP] = getSigmaPointApp_status_mod(varargin)

% Modified version of the getSigmaPointApp.m function of the
% SPToolbox.

nonfun = varargin{1};
xi  = varargin{2};
estruct = varargin{3};
op_SP = varargin{4};
if op_SP.nderiv
    compute_derivative = 1;
else
    compute_derivative = 0;
end

%% Initialization
% Dimensions

beta = estruct.beta(xi);
delta = estruct.delta(xi);

[D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,op_SP.type_D);
n_b = size(D,1);

%% Parameters of sigma point approach
L = n_b;
status = 1;
% Recommendations in:
%       R. van der Merwe, Sigma-point Kalman filters for probabilistic
%       inference in dynamic state-space models, Ph.D. thesis, 2004.
%
% 1) kap >= 0 to guarantee positive semi-definiteness of the covariance
%    matrix. The specific value of kap is not critical though, so a 
%    good default choice is kap= 0.
% 2) Choose 0 <= alp <= 1. 
%    alp controls the "size" of the sigma-point distribution and should
%    ideally be a small number to avoid sampling non-local effects when
%    the nonlinearities are strong.
% 3) Choose bet >= 1. bet is a non-negative weighting term which can be
%    used to incorporate knowledge of the higher order moments of the
%    distribution. For a Gaussian prior the optimal choice is bet = 2.

alp = 0.7;
bet = 2;
kap = 0;
lam = alp^2*(L+kap) - L;

% Weights
w0_m = lam/(L+lam);
wi_m = 1/(2*(L+lam));
w_m = [w0_m;wi_m*ones(2*L,1)];

w0_c = lam/(L+lam)+(1-alp^2+bet);
wi_c = 1/(2*(L+lam));
w_c = [w0_c;wi_c*ones(2*L,1)];

if compute_derivative == 1
    n_xi = size(estruct.ddeltadxi(xi),2);
end

%% Sigma points
% Matrix root: D = S*S'
if compute_derivative == 0
    [S] = chol(D,'lower');
else
    [S,dSddelta] = chol_w_diff(D,dDddelta);
    dSdxi = permute(sum(bsxfun(@times,dSddelta,permute(estruct.ddeltadxi(xi),[3,4,1,2])),3),[1,2,4,3]);
end

% Sigma points
SP.B_SP = [zeros(L,1),sqrt(L+lam)*S,-sqrt(L+lam)*S];
if compute_derivative == 1
    SP.dB_SPdxi = permute([zeros(L,1,n_xi),sqrt(L+lam)*dSdxi,-sqrt(L+lam)*dSdxi],[1,3,2]);
end

%% Propagation of sigma points
% Loop: Sigma points
for i = 1:(2*L+1)
    if compute_derivative == 0
        [status,~,~,Y(:,:,i)] = nonfun(estruct.phi(beta,SP.B_SP(:,i)));
    else
        [status,~,~,Y(:,:,i),~,dYdphi(:,:,:,i)] = nonfun(estruct.phi(beta,SP.B_SP(:,i)));
        dphidxi(:,:,i) = estruct.dphidbeta(beta,SP.B_SP(:,i))*estruct.dbetadxi(xi) + estruct.dphidb(beta,SP.B_SP(:,i))*SP.dB_SPdxi(:,:,i);
        dYdxi(:,:,:,i) = permute(sum(bsxfun(@times,dYdphi(:,:,:,i),permute(dphidxi(:,:,i),[4,3,1,2])),3),[1,2,4,3]);
    end
    if status < 0 
       return;
    end
end

SP.Y = Y;
if(any(isnan(Y)))
    error('Failed to successfully integrate system at all SigmaPoints')
end

[n_t,n_y,~] = size(Y);


%% Evaluation of mean, covariance and cross-covariance
% Mean
if(any([op_SP.req(1),op_SP.req(1),op_SP.req(4),op_SP.req(5)]))
    SP.my = sum(bsxfun(@times,permute(w_m,[3,2,1])  ,Y)    ,3);
    DeltaY = bsxfun(@minus,Y,SP.my);
    if compute_derivative == 1
        SP.dmydxi = sum(bsxfun(@times,permute(w_m,[4,3,2,1]),dYdxi),4);
        DeltadYdxi = bsxfun(@minus,dYdxi,SP.dmydxi);
    end
end

% Covariance
% for measurement noise we ignore the random effects at this point
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
if compute_derivative == 1
    dsigmadphi = estruct.dsigma_noisedphi(estruct.phi(beta,SP.B_SP(:,1)));
    ndim_dsigmadphi = ndims(dsigmadphi);
    dsigmadxi = permute(sum(bsxfun(@times,dsigmadphi,permute(dphidxi(:,:,1),[3:(3+ndim_dsigmadphi-2),1,2])),ndim_dsigmadphi),[1:ndim_dsigmadphi-1,ndim_dsigmadphi+1,ndim_dsigmadphi]);
    dC_techdxi = zeros(n_t,n_y,n_y,size(dsigmadxi,ndim_dsigmadphi));
    
    if(size(dsigmadxi,1) == n_t)
        if(size(dsigmadxi,2) == 1)
            dC_techdxi = bsxfun(@times,repmat(permute(dsigmadxi,[1,2,4,3]),[1,n_y,n_y,1]),permute(eye(n_y),[3,1,2]));
        elseif(size(dsigmadxi,2) == n_y)
            dC_techdxi = bsxfun(@times,repmat(permute(dsigmadxi,[1,2,4,3]),[1,1,n_y,1]),permute(eye(n_y),[3,1,2]));
        end
    elseif(size(dsigmadxi,2) == n_y)
        if(size(dsigmadxi,1) == 1)
            dC_techdxi = bsxfun(@times,repmat(permute(dsigmadxi,[1,2,4,3]),[n_t,1,n_y,1]),permute(eye(n_y),[3,1,2]));
        end
    elseif(and(size(dsigmadxi,1)==1,size(dsigmadxi,2)==1))
        dC_techdxi = bsxfun(@times,repmat(permute(dsigmadxi,[1,2,4,3]),[n_t,n_y,n_y,1]),permute(eye(n_y),[3,1,2]));
    end
end
end


if(op_SP.req(2))
    SP.Cy = sum(bsxfun(@times,permute(w_c,[4,3,2,1]),...
            bsxfun(@times,permute(DeltaY,[1,2,4,3]),permute(DeltaY,[1,4,2,3]))),4) ...
            + C_tech;
    
    if compute_derivative == 1
        dCj =  sum(bsxfun(@times,permute(w_c,[5,4,3,2,1]),...
            bsxfun(@times,permute(DeltaY,[1,2,4,5,3]),permute(DeltadYdxi,[1,5,2,3,4]))),5);
        
        SP.dCydxi = dCj + permute(dCj,[1,3,2,4]) ...
            + dC_techdxi;
        
    end
end

if(op_SP.req(3))
    % Cross-covariance
    SP.Cxy = sum(bsxfun(@times,permute(w_c,[4,3,2,1]),...
        bsxfun(@times,permute(SP.B_SP,[3,1,4,2]),permute(DeltaY,[1,4,2,3]))),4);
    
    if compute_derivative == 1
        SP.dCxydxi = sum(bsxfun(@times,permute(w_c,[5,4,3,2,1]),...
            bsxfun(@times,permute(SP.B_SP,[3,1,4,5,2]),permute(DeltadYdxi,[1,5,2,3,4]))),5) ...
            + sum(bsxfun(@times,permute(w_c,[5,4,3,2,1]),...
                bsxfun(@times,permute(SP.dB_SPdxi,[4,1,5,2,3]),permute(DeltaY,[1,4,2,5,3]))),5);
    end
end

% Full state-time covariance
if(op_SP.req(4))
    SP.mz = SP.my(:);
    if compute_derivative == 1
        SP.dmzdxi = reshape(SP.dmydxi,[size(SP.dmydxi,1)*size(SP.dmydxi,2),size(SP.dmydxi,3)]);
    end
end

if(op_SP.req(5))
    DeltaZ = reshape(DeltaY,[size(DeltaY,1)*size(DeltaY,2),size(DeltaY,3)]);
    SP.Cz = sum(bsxfun(@times,permute(w_c,[3,2,1]),bsxfun(@times,permute(DeltaZ,[1,3,2]),permute(DeltaZ,[3,1,2]))),3) + Cz_tech;
    
    if compute_derivative == 1
        DeltadZdxi = reshape(DeltadYdxi,[size(DeltadYdxi,1)*size(DeltadYdxi,2),size(DeltadYdxi,3),size(DeltadYdxi,4)]);
        dCzdxi = sum(bsxfun(@times,permute(w_c,[4,3,2,1]),bsxfun(@times,permute(DeltaZ,[1,4,3,2]),permute(DeltadZdxi,[4,1,2,3]))),4);
        SP.dCzdxi = dCzdxi + permute(dCzdxi,[2,1,3]);
    end
end

