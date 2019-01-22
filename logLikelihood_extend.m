function varargout = logLikelihood_extend(xi,M,D,varargin)
% This function evaluates the likelihood function for a given model, data
% and parameter vector.
%
% USAGE:
% [...] = logLikelihood(xi,M,D,options,conditions,I) \n
% [...] = logLikelihood(xi,M,D,options,conditions) \n
% [...] = logLikelihood(xi,M,D,options) \n
% [logL] = logLikelihood(...) \n
% [logL, dlogL] = logLikelihood(...) \n
%
% Parameters:
%   xi: parameter values
%   M: model struct
%   D: data struct
%   varargin:
%   options: struct
%   conditions: generated by function collectConditions.m
%   I: indices for which of the data the likelihood function should be evaluated
%
% Return values:
%   logL: log-likelihood value
%   dlogL: gradient of log-likelihood function
%
% Required fields of M:
%   n_subpop: number of subpopulations
%   model: simulation file with input (T,theta,u) (e.g., generated by
%                amiwrap),
%                the first output needs to be the status of the
%                simulation,
%                the 4th the simulation output and
%                the 6th the sensitivities (n_t x n_obs x n_theta)
%   mean_ind: indices of output for mean
%   var_ind: indices of output for variances (empty if using RREs)
%   theta: parameters needed for simulation dependend on xi and u
%       the following fields of M are generated by generate_ODEMM
%   distribution{s,e}: distribution assumption\n
%       = ''norm'': normal distribution assumption\n
%       = ''logn_median'': log-normal distribution assumption, mean of simulation linked to median of distribution\n
%       = ''logn_mean'': log-normal distribution assumption, mean of
%       simulation linked to mean of distribution\n\n
%   The following fields are automatically added by generateODEMM.m
%   dthetadxi: gradient of theta
%   mu{s,e}: specification of mixture parameter mu for subpopulation s and
%                experiment e
%   dmudxi{s,e}: gradient of mu
%   sigma{s,e}: specification of mixture parameter \f$\sigma\f$ (M.Sigma in
%               multivariate case (covariance matrix))
%   dsigmadxi{s,e}: gradient of sigma (M.dSigmadxi in multivariate case)
%   w{s,e}: specification of weights \f$w_s\f$
%   dwdxi{s,e}: gradient of weights
%   scaling{r,e}: scaling parameter of replicate r in experiment e
%   dscalingdxi{r,e}: gradient of scaling
%   offset{r,e}: offset
%   doffsetdxi{r,e}: gradient of offset
%
% Required fields of D:
%   n_dim: dimension of the measurements
%   t: 1 x n_t vector of timepoints
%   u: n_maxu x n_u vector of inputs with n_maxu: maximal number
%                       of inputs simulatenously used
%   y: n_u x n_t x n_cells x n_dim data matrix (only needed if replicates are
%           merged and already scaled), dim is the dimension of the
%           measurement
%   c: n_subpop x (n_u + n_differences) corresponding condition (automatically added by calling collectCondition.m)
%   replicate(r).y:  n_u x n_t x n_cells x n_dim data matrix of replicate r in
%                           experiment e (only needed if individual replicates should be
%                       fitted)
%
% Optional fields of options:
%    use_robust: robust calculation of mixture probability\n
%            = true: uses reformulation (default)\n
%            = false: classical calculation (not recommended)
%    simulate_musigma: true if simulation directly provides ...
%    negLogLikelihood: true if negativev log-likelihood required
%    replicates: true if replicates are fitted individually

%% Set default options
options.use_robust = true;
options.replicates = false;
%options.logPosterior = false;
options.negLogLikelihood = false;
options.simulate_musigma = false;
%% Input assignment
if nargin >= 4
    options = setdefault(varargin{1},options);
end
if nargin >= 5
    conditions = varargin{2};
end
if nargin >= 6
    I = varargin{3};
else
    I = 1:length(D);
end

for e = I
    if options.replicates
        replicates{e} = 1:length(D(e).replicate); % consider replicates individually
    else
        replicates{e} = 1; % consider scaled and merged replicates
    end
end

if ~isfield(M,'w_ind')
    for s = 1:M.n_subpop
        for e = I
            M.w_ind{s,e} = [];
        end
    end
end

%% collect all different conditions for the simulations and simulate them
if nargin < 5 || isempty(conditions)
    [conditions,D] = collectConditions(D,M);
end

for c = 1:length(conditions)
    if nargout >=2
        try
            [status,~,~,X_c{c},~,dXdtheta_c{c}] = M.model(conditions(c).time,...
                M.theta(xi,conditions(c).input),conditions(c).input);
            dXdtheta_c{c} = permute(dXdtheta_c{c},[2,3,1]);
        catch e
            disp(e.message)
            status = -1;
        end
    else
        try
            [status,~,~,X_c{c}] = M.model(conditions(c).time,M.theta(xi,conditions(c).input),conditions(c).input);
        catch e
            disp(e.message)
            status = -1;
        end
    end
    if status < 0
        if options.negLogLikelihood
            varargout{1} =  Inf;
            if nargout >=2
                varargout{2} =  Inf;
            end
        else
            varargout{1} =  -Inf;
            if nargout >=2
                varargout{2} =  -Inf;
            end
        end
        return;
    end
end

%% Evaluation of likelihood function
logL = 0;
dlogL = zeros(length(xi),1);

for e = I % Loop: Experimental conditions
    for d = 1:size(D(e).u,2)
        for r = replicates{e}
            %% get parameters for mixture distribution
            for s = 1:M.n_subpop
                u_dse = [D(e).u(:,d);M.u{s,e}];
                dthetadxi{s} = M.dthetadxi(xi,u_dse);
                t_ind = find(conditions(D(e).c(s,d)).time==D(e).t);
                clear X dXdtheta
                Z = X_c{D(e).c(s,d)}(t_ind,[M.mean_ind{s,e},M.var_ind{s,e},M.w_ind{s,e}]);
                if nargout >= 2
                    dZdtheta = dXdtheta_c{D(e).c(s,d)}([M.mean_ind{s,e},M.var_ind{s,e},M.w_ind{s,e}],:,t_ind);
                end
                if options.simulate_musigma
                    if nargout<2
                        Z = getLognMeanVar(Z,D(e).n_dim);
                    else
                        [Z,dZdtheta] = getLognMeanVar(Z,D(e).n_dim,dZdtheta);
                    end
                end
                % scaling and offset
                X(:,1:D(e).n_dim) = bsxfun(@plus,bsxfun(@times,M.scaling{r,e}(xi,u_dse)',Z(:,1:D(e).n_dim)),...
                    M.offset{r,e}(xi,u_dse)');
                if ~isempty(M.var_ind{s,e})
                    s_temp= M.scaling{r,e}(xi,u_dse);
                    temp = tril(ones(D(e).n_dim,D(e).n_dim));
                    temp(temp==0) = NaN;
                    covscale = (s_temp*s_temp').*temp;
                    covscale = covscale(:);
                    covscale = covscale(~isnan(covscale));
                    for n = 1:(D(e).n_dim*(D(e).n_dim+1))/2
                        X(:,D(e).n_dim+n) = covscale(n)*Z(:,D(e).n_dim+n);
                    end
                end
                switch M.distribution{s,e}
                    case {'logn','logn_median','logn_mean','norm'}
                        if D(e).n_dim == 1
                            sigma{s} = M.sigma{s,e}(D(e).t,X,xi,u_dse);
                            mu{s} = M.mu{s,e}(D(e).t,X,sigma{s},xi,u_dse);
                        else
                            Sigma{s} = M.Sigma{s,e}(D(e).t,X,xi,u_dse);
                            mu{s} = M.mu{s,e}(D(e).t,X,Sigma{s},xi,u_dse);
                        end
                    case 'neg_binomial'
                        rho{s} = M.rho{s,e}(D(e).t,X,xi,u_dse);
                        assert(sum(rho{s}>1)==0,'negative binomial distribution requires variance to be greater than the mean')
                        tau{s} = M.tau{s,e}(D(e).t,X,rho{s},xi,u_dse);
                    case 'students_t'
                        nu{s} = M.nu{s,e}(D(e).t,X,xi,u_dse);
                        Sigma{s} = M.Sigma{s,e}(D(e).t,X,xi,u_dse);
                        mu{s} = M.mu{s,e}(D(e).t,X,Sigma{s},xi,u_dse);
                    case 'skew_t'
                        error('to check')
                        nu{s} = M.nu{s,e}(D(e).t,X,xi,u_dse);
                        mu{s} = M.mu{s,e}(D(e).t,X,xi,u_dse);
                        Sigma{s} = M.Sigma{s,e}(D(e).t,X,xi,u_dse);
                        delta{s} = M.delta{s,e}(D(e).t,X,nu{s},xi,u_dse);
                    case 'skew_norm'
                        delta{s} = M.delta{s,e}(D(e).t,X,xi,u_dse);
                        Sigma{s} = M.Sigma{s,e}(D(e).t,X,delta{s},xi,u_dse);
                        mu{s} = M.mu{s,e}(D(e).t,X,delta{s},xi,u_dse);
                    otherwise
                        error(['Check distribution assumption, provided assumption ''' ...
                            M.distribution{s,e} ''' not covered. Only '...
                            '''neg_binomial'',''students_t'',''logn'',''norm'',''skew_t'',''skew_norm'''])
                end
                w{s} = M.w{s,e}(D(e).t,X,xi,u_dse);
                
                % Derivatives of distribution parameters
                if nargout >= 2
                    dXdxi{s} = zeros(numel([M.mean_ind{s,e},M.var_ind{s,e},M.w_ind{s,e}]),length(xi),length(D(e).t));
                    for k = 1:length(D(e).t)
                        sc = M.scaling{r,e}(xi,u_dse);
                        dsdxi_temp = M.dscalingdxi{r,e}(xi,u_dse);
                        dbdxi_temp = M.doffsetdxi{r,e}(xi,u_dse);
                        for n = 1:D(e).n_dim
                            dXdxi{s}(n,:,k) = sc(n)*dZdtheta(n,:,k)*dthetadxi{s} + ...
                                dsdxi_temp(n,:)*Z(k,n) + dbdxi_temp(n,:);
                        end
                        if ~isempty(M.var_ind{s,e})
                            if D(e).n_dim == 2
                                dcovscaledxi = [2*(dsdxi_temp(1,:)*sc(1));...
                                    (dsdxi_temp(1,:)*sc(2)+dsdxi_temp(2,:)*sc(1));...
                                    2*(dsdxi_temp(2,:)*sc(2))];
                            elseif D(e).n_dim == 1
                                dcovscaledxi = 2*(dsdxi_temp(1,:)*sc(1));
                            else
                                n = 1;
                                for iDim1 = 1:D(e).n_dim
                                    for iDim2 = iDim1:D(e).n_dim
                                        if iDim1 == iDim2
                                            dcovscaledxi(n) = 2*(dsdxi_temp(iDim1,:)*sc(iDim1));
                                        else
                                            dcovscaledxi(n) = dsdxi_temp(iDim1,:)*sc(iDim2)+...
                                                dsdxi_temp(iDim2,:)*sc(iDim1);
                                        end
                                        n=n+1;
                                    end
                                end
                            end
                            for n = 1:(D(e).n_dim*(D(e).n_dim+1))/2
                                dXdxi{s}(D(e).n_dim+n,:,k) = Z(k,D(e).n_dim+n)*(dcovscaledxi(n,:))+...
                                    covscale(n)*dZdtheta(D(e).n_dim+n,:,k)*dthetadxi{s};
                            end
                        end
                    end % time loop
                    switch M.distribution{s,e}
                        case {'logn','logn_median','logn_mean','norm'}
                            if D(e).n_dim == 1
                                dsigmadxi{s} = M.dsigmadxi{s,e}(D(e).t,X,dXdxi{s},xi,u_dse);
                                dmudxi{s} = M.dmudxi{s,e}(D(e).t,X,dXdxi{s},sigma{s},dsigmadxi{s},xi,u_dse);
                            else
                                dSigmadxi{s} = M.dSigmadxi{s,e}(D(e).t,X,dXdxi{s},xi,u_dse);
                                dmudxi{s} = M.dmudxi{s,e}(D(e).t,X,dXdxi{s},Sigma{s},dSigmadxi{s},xi,u_dse);
                            end
                        case 'neg_binomial'
                            drhodxi{s} = M.drhodxi{s,e}(D(e).t,X,dXdxi{s},xi,u_dse);
                            dtaudxi{s} = M.dtaudxi{s,e}(D(e).t,X,dXdxi{s},rho{s},drhodxi{s},xi,u_dse);
                        case 'students_t'
                            dnudxi{s} = M.dnudxi{s,e}(D(e).t,X,dXdxi{s},xi,u_dse);
                            dmudxi{s} = M.dmudxi{s,e}(D(e).t,X,dXdxi{s},xi,u_dse);
                            dSigmadxi{s} = M.dSigmadxi{s,e}(D(e).t,X,dXdxi{s},xi,u_dse);
                        case 'skew_t'
                            error('')
                        case 'skew_norm'
                            ddeltadxi{s} = M.ddeltadxi{s,e}(D(e).t,X,dXdxi{s},xi,u_dse);
                            dSigmadxi{s} = M.dSigmadxi{s,e}(D(e).t,X,dXdxi{s},delta{s},ddeltadxi{s},xi,u_dse);
                            dmudxi{s} = M.dmudxi{s,e}(D(e).t,X,dXdxi{s},delta{s},ddeltadxi{s},xi,u_dse);
                    end
                    dwdxi{s} = M.dwdxi{s,e}(D(e).t,X,dXdxi{s},xi,u_dse);
                end % gradient
            end % subpopulation
            
            % Loop over the time points and
            for k = 1:length(D(e).t)
                % get data
                if options.replicates
                    y = squeeze(D(e).replicate(r).y(d,k,:,:));
                else
                    y = squeeze(D(e).y(d,k,:,:));
                end
                y = y((sum(~isnan(y),2) == size(y,2)),:);
                
                % initialize
                if options.use_robust
                    q = zeros(length(y),M.n_subpop);
                    w_s = [];
                    if nargout >= 2
                        H = zeros(length(y),length(xi),M.n_subpop);
                    end
                else
                    p = zeros(length(y),1);
                    dpdxi = zeros(length(y),length(xi));
                end
                %% evaluate likelihood function components
                for s = 1:M.n_subpop
                    if options.use_robust
                        switch M.distribution{s,e}
                            case {'logn','logn_median','logn_mean'}
                                if D(e).n_dim == 1
                                    q(:,s) = logoflognpdf(y,mu{s}(k), sigma{s}(k));
                                    if nargout >= 2
                                        H(:,:,s) = bsxfun(@plus,dwdxi{s}(k,:),w{s}(k)/sigma{s}(k)*...
                                            ((log(y)-mu{s}(k))/sigma{s}(k)*dmudxi{s}(k,:)+...
                                            (((log(y)-mu{s}(k))/sigma{s}(k)).^2-1)*dsigmadxi{s}(k,:)));
                                    end
                                else % multivariate
                                    q(:,s) = bsxfun(@minus,logofmvnpdf(log(y),mu{s}(k,:),permute(Sigma{s}(k,:,:),[2,3,1])),sum(log(y),2));
                                    if nargout >= 2
                                        if rcond(permute(Sigma{s}(k,:,:),[2 3 1])) < 1e-10
                                            error('Sigma bad scaled')
                                        end
                                        SigmaIn = inv(permute(Sigma{s}(k,:,:),[2 3 1]));
                                        for n_xi = 1:length(xi)
                                            dSigmaIndxi = -SigmaIn*permute(dSigmadxi{s}(k,n_xi,:,:),[3,4,1,2])*SigmaIn;
                                            H(:,n_xi,s) = bsxfun(@plus,dwdxi{s}(k,n_xi), w{s}(k).*(-0.5)*...
                                                (repmat(sum(sum((SigmaIn.').*permute(dSigmadxi{s}(k,n_xi,:,:),[3,4,1,2]))),length(y),1)...
                                                +  bsxfun(@minus,mu{s}(k,:),log(y))*SigmaIn*permute(dmudxi{s}(k,n_xi,:),[3,1,2])...
                                                + (permute(dmudxi{s}(k,n_xi,:),[1,3,2])*SigmaIn*(bsxfun(@minus,mu{s}(k,:),log(y)))')'...
                                                + sum((bsxfun(@minus,mu{s}(k,:),log(y))*dSigmaIndxi).*bsxfun(@minus,mu{s}(k,:),log(y)),2)));
                                        end
                                    end
                                end
                            case 'norm'
                                if D(e).n_dim == 1
                                    q(:,s) = logofnormpdf(y,mu{s}(k), sigma{s}(k));
                                    if nargout >= 2
                                        H(:,:,s) = bsxfun(@plus,dwdxi{s}(k,:),w{s}(k)/sigma{s}(k)*...
                                            ((y-mu{s}(k))/sigma{s}(k)*dmudxi{s}(k,:)+...
                                            (((y-mu{s}(k))/sigma{s}(k)).^2-1)*dsigmadxi{s}(k,:)));
                                    end
                                else % multivariate
                                    q(:,s) = logofmvnpdf(y,mu{s}(k,:),permute(Sigma{s}(k,:,:),[2,3,1]));
                                    if nargout >= 2
                                        SigmaIn = inv(permute(Sigma{s}(k,:,:),[2 3 1]));
                                        for n_xi = 1:length(xi)
                                            dSigmaIndxi = -SigmaIn*permute(dSigmadxi{s}(k,n_xi,:,:),[3,4,1,2])*SigmaIn;
                                            
                                            H(:,n_xi,s) = bsxfun(@plus,dwdxi{s}(k,n_xi), w{s}(k).*(-0.5)*...
                                                (repmat(sum(sum((SigmaIn.').*permute(dSigmadxi{s}(k,n_xi,:,:),[3,4,1,2]))),length(y),1)...
                                                +  bsxfun(@minus,mu{s}(k,:),y)*SigmaIn*permute(dmudxi{s}(k,n_xi,:),[3,1,2])...
                                                + (permute(dmudxi{s}(k,n_xi,:),[1,3,2])*SigmaIn*(bsxfun(@minus,mu{s}(k,:),y))')'...
                                                + sum((bsxfun(@minus,mu{s}(k,:),y)*dSigmaIndxi).*bsxfun(@minus,mu{s}(k,:),y),2)));
                                            
                                        end % xi
                                    end % gradient
                                end % dimension
                            case 'neg_binomial'
                                if nargout<2
                                    q(:,s) = logofnbinpdf(y,tau{s}(k),rho{s}(k));
                                else
                                    [q(:,s),dqdxi] = logofnbinpdf(y,tau{s}(k),rho{s}(k),dtaudxi{s}(k,:),drhodxi{s}(k,:));
                                    H(:,:,s) = bsxfun(@plus,dwdxi{s}(k,:),w{s}(k)*dqdxi);
                                end
                            case 'students_t'
                                if nargout<2
                                    q(:,s) = logofmvtpdf(y,mu{s}(k,:),permute(Sigma{s}(k,:,:),[2,3,1]),nu{s}(k));
                                else
                                    [q(:,s),dqdxi] = logofmvtpdf(y,mu{s}(k,:),permute(Sigma{s}(k,:,:),[2,3,1]),nu{s}(k),...
                                        permute(dmudxi{s}(k,:,:),[3,2,1]),permute(dSigmadxi{s}(k,:,:,:),[3,4,1,2]),dnudxi{s}(k,:));
                                    H(:,:,s) = bsxfun(@plus,dwdxi{s}(k,:),w{s}(k)*dqdxi');
                                end
                            case 'skew_norm'
                                if nargout<2
                                    q(:,s) = logofskewnormpdf(y,mu{s}(k,:),...
                                        permute(Sigma{s}(k,:,:),[2,3,1]),delta{s});
                                else
                                    [q(:,s),dqdxi] = logofskewnormpdf(y,mu{s}(k,:),permute(Sigma{s}(k,:,:),[2,3,1]),delta{s},...
                                        permute(dmudxi{s}(k,:,:),[3,2,1]),permute(dSigmadxi{s}(k,:,:,:),[3,4,1,2]),ddeltadxi{s});
                                    H(:,:,s) = bsxfun(@plus,dwdxi{s}(k,:),w{s}(k)*dqdxi');
                                end
                        end % distribution
                        w_s = [w_s,w{s}(k)];
                    else % not robust, not recommended
                        switch M.distribution{s,e}
                            case {'logn','logn_median','logn_mean'}
                                if D(e).n_dim == 1
                                    p_s = pdf('logn',y,mu{s}(k),sigma{s}(k));
                                    p = p + w{s}(k)*p_s;
                                    if nargout >= 2
                                        dpdxi = dpdxi + p_s*dwdxi{s}(k,:) + ...
                                            w{s}(k)/sigma{s}(k)*bsxfun(@times,p_s,...
                                            ((log(y)-mu{s}(k))/sigma{s}(k)*dmudxi{s}(k,:)+...
                                            (((log(y)-mu{s}(k))/sigma{s}(k)).^2-1)*dsigmadxi{s}(k,:)));
                                        
                                    end
                                else % multivariate
                                    p_s = bsxfun(@rdivide,mvnpdf(log(y),mu{s}(k,:),permute(Sigma{s}(k,:,:),[2,3,1])),prod(y,2));
                                    p = p + w{s}(k)*p_s;
                                    if nargout >= 2
                                        if rcond(permute(Sigma{s}(k,:,:),[2 3 1])) < 1e-10
                                            error('Sigma bad scaled')
                                        end
                                        SigmaIn = inv(permute(Sigma{s}(k,:,:),[2 3 1]));
                                        permute(Sigma{s}(k,:,:),[2 3 1])
                                        for n_xi = 1:length(xi)
                                            dSigmaIndxi = -SigmaIn*permute(dSigmadxi{s}(k,n_xi,:,:),[3,4,1,2])*SigmaIn;
                                            dpdxi(:,n_xi) = dpdxi(:,n_xi) + p_s.*dwdxi{s}(n_xi) +...
                                                + w{s}(k).*(-0.5).*p_s.*(repmat(sum(sum((SigmaIn.').*permute(dSigmadxi{s}(k,n_xi,:,:),[3,4,1,2]))),length(y),1)...
                                                + bsxfun(@minus,mu{s}(k,:),log(y))*SigmaIn*permute(dmudxi{s}(k,n_xi,:),[3,1,2])...
                                                + (permute(dmudxi{s}(k,n_xi,:),[1,3,2])*SigmaIn*(bsxfun(@minus,mu{s}(k,:),log(y)))')'...
                                                + sum((bsxfun(@minus,mu{s}(k,:),log(y))*dSigmaIndxi).*bsxfun(@minus,mu{s}(k,:),log(y)),2));
                                        end
                                    end
                                end
                            case 'norm'
                                if D(e).n_dim == 1
                                    p_s = pdf('norm',y,mu{s}(k),sigma{s}(k));
                                    p = p + w{s}(k)*p_s;
                                    if nargout >= 2
                                        dpdxi = dpdxi + p_s*dwdxi{s}(k,:) + ...
                                            w{s}(k)/sigma{s}(k)*bsxfun(@times,p_s,...
                                            ((y-mu{s}(k))/sigma{s}(k)*dmudxi{s}(k,:)+...
                                            (((y-mu{s}(k))/sigma{s}(k)).^2-1)*dsigmadxi{s}(k,:)));
                                    end
                                else % multivariate
                                    p_s = mvnpdf(y,mu{s}(k,:),permute(Sigma{s}(k,:,:),[2,3,1]));
                                    p = p + w{s}(k)*p_s;
                                    if nargout >= 2
                                        SigmaIn = inv(permute(Sigma{s}(k,:,:),[2 3 1]));
                                        for n_xi = 1:length(xi)
                                            dSigmaIndxi = -SigmaIn*permute(dSigmadxi{s}(k,n_xi,:,:),[3,4,1,2])*SigmaIn;
                                            dpdxi(:,n_xi) = dpdxi(:,n_xi) + p_s.*dwdxi{s}(k,n_xi) +...
                                                + w{s}(k).*(-0.5).*p_s.*(repmat(sum(sum((SigmaIn.').*permute(dSigmadxi{s}(k,n_xi,:,:),[3,4,1,2]))),length(y),1)...
                                                + bsxfun(@minus,mu{s}(k,:),y)*SigmaIn*permute(dmudxi{s}(k,n_xi,:),[3,1,2])...
                                                + (permute(dmudxi{s}(k,n_xi,:),[1,3,2])*SigmaIn*(bsxfun(@minus,mu{s}(k,:),y))')'...
                                                + sum((bsxfun(@minus,mu{s}(k,:),y)*dSigmaIndxi).*bsxfun(@minus,mu{s}(k,:),y),2));
                                            
                                        end % xi
                                    end % gradient
                                end % dimension
                            case {'neg_binomial'}
                                p_s = nbinpdf(y,tau{s}(k),rho{s}(k));
                                p = p + w{s}(k)*p_s;
                                if nargout >= 2
                                    dpdxi = dpdxi + p_s.*dwdxi{s}(k,:) +...
                                        + w{s}(k)*(bsxfun(@times,psi(y+tau{s}(k)),dtaudxi{s}(k,:))-...
                                        psi(tau{s}(k))*dtaudxi{s}(k) + ...
                                        dtaudxi{s}(k,:)*log(1-rho{s}(k)) -...
                                        y/(1-rho{s}(k,:))*drhodxi{s}(k,:)+...
                                        tau{s}(k)./rho{s}(k).*drhodxi{s}(k,:));
                                end
                            case 'student_t'
                                error('non-robust version for student t not supported')
                        end % distribution
                    end % robust
                end % subpopulation loop
                
                %% Evaluation of mixture likelihood for this time point and dose
                if options.use_robust
                    if nargout >= 2
                        [logp,dlogpdxi]= computeMixtureProbability(w_s,q,H) ;
                        if size(dlogpdxi,1)>1
                            dlogL = dlogL + sum(dlogpdxi)';
                        else
                            dlogL = dlogL + dlogpdxi';
                        end
                    elseif nargout <= 1
                        logp = computeMixtureProbability(w_s,q) ;
                    end
                else
                    logp = log(p);
                    if nargout >= 2
                        dlogL = dlogL + sum(bsxfun(@times,1./p,dpdxi))';
                    end
                end
                logL = logL + sum(logp);
            end % time loop
        end % replicates
    end % dose loop
end % experiment


J = logL;

%% Output assignment
if ~isreal(J)
    error('Likelihood is not real!');
else
    if options.negLogLikelihood
        varargout{1} =  -J;
    else
        varargout{1} =  J;
    end
end
if nargout >=2
    if ~isreal(dlogL)
        error('Gradient is not real');
    elseif sum(isnan(dlogL))>0
        error('Gradient contains NaNs');
    elseif sum(isinf(dlogL))>0
        error('Gradient contains Infs');
    else
        if options.negLogLikelihood
            varargout{2} =  -dlogL;
        else
            varargout{2} =  dlogL;
        end
    end
end

