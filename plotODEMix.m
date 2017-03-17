% function [fh,fhm] = plotODEMix(D,M,xi,I,options,tu_ind)
function varargout = plotODEMix(varargin)
% Routine to plot the ODE-constrained mixture model
%
% USAGE:
% [...] = plotODEMix(D,M,xi) \n
% [...] = plotODEMix(D,M,xi,I,options) \n
% [...] = plotODEMix(D,M,xi,I,options,tu_ind) \n
% [fh] = plotODEMix(...) \n
% [fh,fhm] = plotODEMix(...)
%
% Parameters:
% varargin:
% * D: data struct
% * M: model struct 
% * xi: parameter vector
% * I: (optional) indices for which the data and model should be visualized,
%   the whole data set is visualized if I = []
% * options: plotting options
% * tu_ind{e}: struct of indices for the time points/doses for which the data and model
% should be visualized
%
% Return values:
% fh: struct of function handles for each data set
% fhm: struct of function handles for the plots of the marginals


%% Assign inputs
D = varargin{1};
if nargin >= 3
    M = varargin{2};
    xi = varargin{3};
else
    M = [];
    xi = [];
end
I = [];
if nargin >=4
    I = varargin{4};
end
if isempty(I)
    I = 1:length(D);
end

if nargin >= 6
    tu_ind = varargin{6};
    if isempty(tu_ind)
        for e=I
            if numel(D(e).t) > 1
                tu_ind{e} = [1:numel(D(e).t)];
            else
                tu_ind{e} = 1:size(D(e).u,2);
            end
        end
    end
else
    for e=I
        if numel(D(e).t) > 1
            tu_ind{e} = [1:numel(D(e).t)];
        else
            tu_ind{e} = 1:size(D(e).u,2);
        end
    end
end


% Set defaults
options.type = 'kde';
options.switch_axes = false;
options.hold_on = 0;
% model
options.model.col = 'r';
options.model.lw  = 2;
options.model.points = 100;
options.model.plot = 'continuous'; % 'hist'
options.model.subpopulations = false;
options.model.ls = '-';
% 2D
options.model.levelsets = 5;
options.model.level_linewidth = 1.2;
options.model.colormap = 'autumn';
% data
options.data.kde = false;
options.data.plot = 'filled'; %'empty'
options.data.col  = 'b';
options.data.fill_col = 0.7*[1,1,1];
options.data.lw   = 2;
options.data.bins = 100;
% 2D
options.data.marker = 'k.';
options.data.markersize = 1.5;

options.marginals = false;
options.replicates = false;
options.x_scale = 'lin';
options.simulate_musigma = false;
options.sameplot = false;
options.subplot_lin = false;
options.plainstyle = false;
options.legendflag = true;

for e = 1:length(D)
    options.boundaries(e).y_min = []; % 1 x n_meas vector
    options.boundaries(e).y_max = []; % 1 x n_meas vector
    options.ylabel{e} = 'frequency';
    options.z_max{e} = [];
    options.xtick{e} = [];
    options.ytick{e} = [];
end
if nargin >= 5
    options = setdefault(varargin{5},options);
end
if ~iscell(options.model.col)
    temp = options.model.col;
    options.model = rmfield(options.model, 'col');
    for e = I
        options.model.col{e} = temp;
    end
end

if ~iscell(options.data.col)
    temp = options.data.col;
    options.data = rmfield(options.data, 'col');
    for e = I
        options.data.col{e} = temp;
    end
end

if ~iscell(options.data.fill_col)
    temp = options.data.fill_col;
    options.data = rmfield(options.data, 'fill_col');
    for e = I
        options.data.fill_col{e} = temp;
    end
end

if ~iscell(options.model.colormap)
    temp = options.model.colormap;
    options.model = rmfield(options.model, 'colormap');
    for e = I
        options.model.colormap{e} = temp;
    end
end

for e = I
    if options.replicates
        n_replicates{e} = 1:length(D(e).replicate); % consider replicates individually
    else
        n_replicates{e} = 1; % consider scaled and merged replicates
    end
end
if options.hold_on
    fh = varargin{7};
end
%% simulate conditions
if ~isempty(M)
    [conditions,D] = collectConditions(D,M);
    for c = 1:length(conditions)
        [~,~,~,X_c{c}] = M.model(conditions(c).time,M.theta(xi,conditions(c).input),conditions(c).input);
    end
end


%% loop over experimental conditions that need to be plotted
for e = I
    clearvars w mu sigma Sigma nu
    r=1;
    %% check plotting case
    if size(D(e).u,2) == 1 && size(D(e).t,2) >= 1
        % one dose with (one or more) time points
        plotcase = 'one dose more tps';
    elseif size(D(e).u,2) >= 1 && size(D(e).t,2) == 1
        % one time point with more dose
        plotcase = 'more dose one tp';
    end
    %% open figures
    if (~options.sameplot || mod(e,2)) & ~options.hold_on
        
        if D(e).n_dim > 1 && options.marginals
            for n = 1:D(e).n_dim
                fhm{e,n} = figure('name',[D(e).name ...
                    ', marginal for ' D(e).measurand{n} '']);
            end
        end
        fh{e} = figure('name',[D(e).name]);
        
    end
    %% evaluate model
    inds = 0;
    if options.marginals && D(e).n_dim > 1
        inds = [0:D(e).n_dim];
    end
    for ind = 1:numel(inds)
        switch plotcase
            case 'one dose more tps'
                if inds(ind) > 0
                    figure(fhm{e,inds(ind)});
                end
                d=1;
                c=1;
                if ~isempty(M)
                    evalModel(xi,M,D,e,r,d,X_c,options,conditions);
                end
                for k=1:numel(tu_ind{e})
                    [lim,hists,grids]=setYminmaxHists(D,e,tu_ind{e}(d),options,inds(ind));
                    if D(e).n_dim == 2 && ~isempty(M)
                        if inds(ind) == 0
                            if options.data.kde
                                for i = 1:2
                                    subplot(2,numel(tu_ind{e}),c);
                                    if i == 1
                                        c = c+numel(tu_ind{e});
                                    else
                                        c = c-numel(tu_ind{e})+1;
                                    end
                                    evalPdf(M,D,e,d,tu_ind{e}(k),options,...
                                        (options.legendflag & k==numel(tu_ind{e})),...
                                        inds(ind),lim,hists,grids,i==2,i==1);
                                end
                            else
                                sx = round(sqrt(numel(tu_ind{e})));
                                sy = ceil(numel(tu_ind{e})/sx);
                                subplot(sx,sy,k);
                                evalPdf(M,D,e,d,tu_ind{e}(k),options,...
                                    (options.legendflag & k==numel(tu_ind{e})),...
                                    inds(ind),lim,hists,grids,1,0);
                            end
                        else
                            if options.subplot_lin
                                subplot(1,numel(tu_ind{e}),k);
                            else
                                sx = round(sqrt(numel(tu_ind{e})));
                                sy = ceil(numel(tu_ind{e})/sx);
                                subplot(sx,sy,k);
                            end
                            c=c+1;
                            evalPdf(M,D,e,d,tu_ind{e}(k),options,...
                                (options.legendflag & k==numel(tu_ind{e})),...
                                inds(ind),lim,hists,grids,1,1);
                        end
                    elseif D(e).n_dim == 2 && isempty(M)
                        
                        if options.subplot_lin
                            subplot(1,numel(tu_ind{e}),k);
                        else
                            sx = round(sqrt(numel(tu_ind{e})));
                            sy = ceil(numel(tu_ind{e})/sx);
                            subplot(sx,sy,k);
                        end
                        evalPdf(M,D,e,d,tu_ind{e}(k),options,...
                            (options.legendflag & k==numel(tu_ind{e})),...
                            inds(ind),lim,hists,grids,0,1);
                    else
                        if options.subplot_lin
                            subplot(1,numel(tu_ind{e}),k);
                        else
                            sx = round(sqrt(numel(tu_ind{e})));
                            sy = ceil(numel(tu_ind{e})/sx);
                            subplot(sx,sy,k);
                        end
                        evalPdf(M,D,e,d,tu_ind{e}(k),options,...
                            (options.legendflag & k==numel(tu_ind{e})),...
                            inds(ind),lim,hists,grids,1,1);
                        if tu_ind{e}(k)~=1
                            if options.plainstyle
                                set(gca,'xtick','');
                                axis off
                            end
                        end
                        if options.plainstyle
                            set(gca,'ytick','','XMinorTick','off');
                        end
                        
                        set(gca,'ticklength',2*get(gca,'ticklength'))
                        box off
                        if options.subplot_lin
                            view(90,-90);
                        end
                    end
                end
            case 'more dose one tp'
                if inds(ind) > 0
                    figure(fhm{e,inds(ind)});
                else
                    figure(fh{e});
                end
                k=1;
                c=1;
                for d=1:numel(tu_ind{e})
                    if ~isempty(M)
                        evalModel(xi,M,D,e,r,tu_ind{e}(d),X_c,options,conditions);
                    end
                    [lim,hists,grids]=setYminmaxHists(D,e,tu_ind{e}(d),options,inds(ind));

                    if D(e).n_dim == 2 && ~isempty(M)
                        if inds(ind) == 0
                            if options.data.kde
                                for i = 1:2
                                    subplot(2,numel(tu_ind{e}),c);
                                    if i == 1
                                        c = c+numel(tu_ind{e});
                                    else
                                        c = c-numel(tu_ind{e})+1;
                                    end
                                    evalPdf(M,D,e,tu_ind{e}(d),k,options,...
                                        (options.legendflag & d==numel(tu_ind{e})),...
                                        inds(ind),lim,hists,grids,i==2,i==1);
                                end
                            else
                                sx = round(sqrt(numel(tu_ind{e})));
                                sy = ceil(numel(tu_ind{e})/sx);
                                subplot(sx,sy,d);
                                evalPdf(M,D,e,tu_ind{e}(d),k,options,...
                                    (options.legendflag & d==numel(tu_ind{e})),...
                                    inds(ind),lim,hists,grids,1,0);
                            end 
                        else
                            if options.subplot_lin
                                subplot(1,numel(tu_ind{e}),d);
                            else
                                sx = round(sqrt(numel(tu_ind{e})));
                                sy = ceil(numel(tu_ind{e})/sx);
                                subplot(sx,sy,d);
                            end
                            c=c+1;
                            evalPdf(M,D,e,tu_ind{e}(d),k,options,0,inds(ind),lim,hists,grids,1,1);
                        end
                    elseif D(e).n_dim == 2 && isempty(M)
                        if options.subplot_lin
                            subplot(1,numel(tu_ind{e}),d);
                        else
                            sx = round(sqrt(numel(tu_ind{e})));
                            sy = ceil(numel(tu_ind{e})/sx);
                            subplot(sx,sy,d);
                        end
                        
                        evalPdf(M,D,e,tu_ind{e}(d),k,options,0,inds(ind),lim,hists,grids,0,1);
                    else
                        if options.subplot_lin
                            subplot(1,numel(tu_ind{e}),d);
                        else
                            sx = round(sqrt(numel(tu_ind{e})));
                            sy = ceil(numel(tu_ind{e})/sx);
                            subplot(sx,sy,d);
                        end
                        
                        evalPdf(M,D,e,tu_ind{e}(d),k,options,0,inds(ind),lim,hists,grids,1,1);
                        
                        if tu_ind{e}(d)~=1
                            if options.plainstyle
                                set(gca,'xtick','');
                                axis off
                            end
                        end
                        if options.plainstyle
                            set(gca,'ytick','','XMinorTick','off');
                        end
                        set(gca,'ticklength',2*get(gca,'ticklength'))
                        box off
                        if options.subplot_lin
                            view(90,-90);
                        end
                    end
                end
            otherwise
                error('Plotcase not clear!')
        end
    end
end
if ~options.hold_on
    varargout{1} = fh;
    if nargout >= 2
        varargout{2} = fhm;
    end
end
end

function str_dose = getStrDose(D,e,d)
str_dose = ['dose = '];
if size(D(e).u,1) > 1
    str_dose = [str_dose '('];
end
for i = 1:size(D(e).u,1)
    str_dose = [str_dose '' num2str(D(e).u(i,d)) ''];
    if i < size(D(e).u,1)
        str_dose = [str_dose ','];
    end
end
if size(D(e).u,1) > 1
    str_dose = [str_dose ')'];
end
end

function [] = evalModel(xi,M,D,e,r,d,X_c,options,conditions)
global w mu sigma Sigma nu
for s = 1:M.n_subpop
    u_dse = [D(e).u(:,d);M.u{s,e}];
    t_ind = find(conditions(D(e).c(s,d)).time==D(e).t);
    clear X dXdtheta
    Z = X_c{D(e).c(s,d)}(t_ind,[M.mean_ind{s,e},M.var_ind{s,e},M.w_ind{s,e}]);
    if options.simulate_musigma
        Z = getLognMeanVar(Z,D(e).n_dim);
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
    if D(e).n_dim == 1
        switch M.distribution{s,e}
            case {'norm', 'logn', 'logn_median', 'logn_mean'}
                sigma{s} = M.sigma{s,e}(D(e).t,X,xi,u_dse);
                mu{s} = M.mu{s,e}(D(e).t,X,sigma{s},xi,u_dse);
            case 't'
                nu{s} = M.nu{s,e}(D(e).t,X,xi,u_dse);
                mu{s} = M.mu{s,e}(D(e).t,X,xi,u_dse);
                sigma{s} = M.sigma{s,e}(D(e).t,X,xi,nu{s},u_dse);
            otherwise
                error('Invalid distribution assumption')
        end
    else
        switch M.distribution{s,e}
            case {'norm', 'logn', 'logn_median', 'logn_mean'}
                Sigma{s} = M.Sigma{s,e}(D(e).t,X,xi,u_dse);
                mu{s} = M.mu{s,e}(D(e).t,X,Sigma{s},xi,u_dse);
            case 't'
                error('todo multidim tdist');
            otherwise
                error('Invalid distribution assumption')
        end
    end
    w{s} = M.w{s,e}(D(e).t,X,xi,u_dse);
end
end

function [lim,hists,grids] = setYminmaxHists(D,e,d,options,ind)
if D(e).n_dim == 1 || ind > 0
    if isempty(options.boundaries(e).y_min) || isempty(options.boundaries(e).y_max)
        y_min{e} =  inf;
        y_max{e} = -inf;
        if (~options.replicates && ~isempty(D(e).y)) || ...
                (options.replicates && ~isempty(D(e).replicate(r).y))
            for k = 1:length(D(e).t)
                if options.replicates
                    if ind > 0
                        y = squeeze(D(e).replicate(r).y(d,k,:,ind));
                    else
                        y = squeeze(D(e).replicate(r).y(d,k,:));
                    end
                else
                    if ind > 0
                        y = squeeze(D(e).y(d,k,:,ind));
                    else
                        y = squeeze(D(e).y(d,k,:));
                    end
                end
                y = y(~isnan(y));
                if ~isempty(y)
                    y_min{e} = min(y_min{e},min(y));
                    y_max{e} = max(y_max{e},max(y));
                end
            end
        else
            error('specify options.boundaries(e).y_min and options.boundaries(e).y_max')
        end
    elseif ind > 0
        y_min{e} = options.boundaries(e).y_min(ind);
        y_max{e} = options.boundaries(e).y_max(ind);
    else
        y_min{e} = options.boundaries(e).y_min;
        y_max{e} = options.boundaries(e).y_max;
    end
    switch options.x_scale
        case 'lin'
            y_hist = linspace(y_min{e},y_max{e},options.data.bins+1)';
            d_y_hist = y_hist(2)-y_hist(1);
            y_grid = linspace(y_min{e},y_max{e},options.model.points+1)';
            d_y_grid = y_grid(2)-y_grid(1);
        case 'log'
            z_hist = linspace(log10(y_min{e}),log10(y_max{e}),options.data.bins+1)';
            d_y_hist = (z_hist(2)-z_hist(1));
            y_hist = 10.^z_hist;
            z_grid = linspace(log10(y_min{e}),log10(y_max{e}),options.model.points+1)';
            d_y_grid = (z_grid(2)-z_grid(1));
            y_grid = 10.^z_grid;
    end
elseif D(e).n_dim == 2 && ind == 0
    if isempty(options.boundaries(e).y_min) || isempty(options.boundaries(e).y_max)
        y_min{e}(1) =  inf;
        y_max{e}(1) = -inf;
        y_min{e}(2) =  inf;
        y_max{e}(2) = -inf;
        for k = 1:length(D(e).t)
            y = squeeze(D(e).y(d,k,:,:));
            y_min{e}(1) = min(y_min{e}(1),min(y(:,1)));
            y_max{e}(1) = max(y_max{e}(1),max(y(:,1)));
            y_min{e}(2) = min(y_min{e}(1),min(y(:,2)));
            y_max{e}(2) = max(y_max{e}(1),max(y(:,2)));
        end
    else
        y_min{e}(1) = options.boundaries(e).y_min(1);
        y_max{e}(1) = options.boundaries(e).y_max(1);
        y_min{e}(2) = options.boundaries(e).y_min(2);
        y_max{e}(2) = options.boundaries(e).y_max(2);
    end % boundaries
    switch options.x_scale
        case 'lin'
            y_hist{1} = linspace(y_min{e}(1),y_max{e}(1),options.data.bins+1)';
            d_y_hist{1} = y_hist{1}(2)-y_hist{1}(1);
            y_grid{1} = linspace(y_min{e}(1),y_max{e}(1),options.model.points+1)';
            d_y_grid{1} = y_grid{1}(2)-y_grid{1}(1);
            y_hist{2} = linspace(y_min{e}(2),y_max{e}(2),options.data.bins+1)';
            d_y_hist{2}  = y_hist{2}(2)-y_hist{2} (1);
            y_grid{2}  = linspace(y_min{e}(2),y_max{e}(2),options.model.points+1)';
            d_y_grid{2}  = y_grid{2}(2)-y_grid{2}(1);
        case 'log'
            z_hist{1} = linspace(log10(y_min{e}(1)),log10(y_max{e}(1)),options.data.bins+1)';
            d_y_hist{1} = (z_hist{1}(2)-z_hist{1}(1));
            y_hist{1} = 10.^z_hist{1};
            z_grid{1} = linspace(log10(y_min{e}(1)),log10(y_max{e}(1)),options.model.points+1)';
            d_y_grid{1} = (z_grid{1}(2)-z_grid{1}(1));
            y_grid{1} = 10.^z_grid{1};
            z_hist{2}  = linspace(log10(y_min{e}(2)),log10(y_max{e}(2)),options.data.bins+1)';
            d_y_hist{2}  = (z_hist{2}(2)-z_hist{2}(1));
            y_hist{2}  = 10.^z_hist{2};
            z_grid{2}  = linspace(log10(y_min{e}(2)),log10(y_max{e}(2)),options.model.points+1)';
            d_y_grid{2}  = (z_grid{2}(2)-z_grid{2}(1));
            y_grid{2}  = 10.^z_grid{2} ;
    end
else
    error('plot only for dimension up to 2');
end
lim.y_min = y_min;
lim.y_max = y_max;
hists.y_hist = y_hist;
hists.d_y_hist = d_y_hist;
grids.y_grid = y_grid;
grids.d_y_grid = d_y_grid;
if strcmp(options.x_scale,'log')
    grids.z_grid = z_grid;
    hists.z_hist = z_hist;
end
end


function [] = evalPdf(M,D,e,d,k,options,legendflag,ind,lim,hists,grids,plotModel,plotData)
global w mu sigma Sigma nu
y_min = lim.y_min;
y_max = lim.y_max;
y_hist = hists.y_hist;
y_grid = grids.y_grid;
d_y_grid = grids.d_y_grid;
d_y_hist = hists.d_y_hist;

if (~options.replicates && ~isempty(D(e).y)) || ...
        (options.replicates && ~isempty(D(e).replicate(r).y))
    if options.replicates
        if ind > 0
            y = squeeze(D(e).replicate(r).y(d,k,:,ind));
        else
            y = squeeze(D(e).replicate(r).y(d,k,:,:));
        end
    else
        if ind > 0
            y = squeeze(D(e).y(d,k,:,ind));
        else
            y = squeeze(D(e).y(d,k,:,:));
        end
        y = y((sum(~isnan(y),2) == size(y,2)),:);
        if D(e).n_dim == 1 || ind > 0
            h = hist(y,y_hist);
            h = (h(1:end-1)'/sum(h(1:end-1)));
        end
    end
end
% Density calculation
if ~isempty(M)
    if (D(e).n_dim == 1) || (ind > 0)
        p = zeros(size(y_grid));
        cp = zeros(size(y_grid));
        for s = 1:M.n_subpop
            switch M.distribution{s,e}
                case {'logn','logn_median','logn_mean'}
                    if ind > 0
                        Sigma_temp = permute(Sigma{s}(k,:,:),[2,3,1]);
                        sigma{s}(k) = sqrt(Sigma_temp(ind,ind));
                    end
                    if ind > 0
                        p = p + w{s}(k)*pdf('logn',y_grid,mu{s}(k,ind),sigma{s}(k));
                        p_s{s} = w{s}(k)*pdf('logn',y_grid,mu{s}(k,ind),sigma{s}(k));
                        cp = cp + w{s}(k)*cdf('logn',y_grid,mu{s}(k,ind),sigma{s}(k));
                    else
                        p = p + w{s}(k)*pdf('logn',y_grid,mu{s}(k),sigma{s}(k));
                        p_s{s} = w{s}(k)*pdf('logn',y_grid,mu{s}(k),sigma{s}(k));
                        cp = cp + w{s}(k)*cdf('logn',y_grid,mu{s}(k),sigma{s}(k));
                    end
                case 'norm'
                    if ind > 0
                        Sigma_temp = permute(Sigma{s}(k,:,:),[2,3,1]);
                        sigma{s}(k) = sqrt(Sigma_temp(n));
                    end
                    p = p + w{s}(k)*pdf('norm',y_grid,mu{s}(k),sigma{s}(k));
                    p_s{s} = w{s}(k)*pdf('norm',y_grid,mu{s}(k),sigma{s}(k));
                    cp = cp + w{s}(k)*cdf('norm',y_grid,mu{s}(k),sigma{s}(k));
                case 't'
                    p = p + w{s}(k)*exp(mylogofmvtpdf(y_grid,mu{s}(k),sigma{s}(k),nu{s}));
                    p_s{s} = w{s}(k)*exp(mylogofmvtpdf(y_grid,mu{s}(k),sigma{s}(k),nu{s}));
            end
        end
    else
        [Y1,Y2] = meshgrid(y_grid{1},y_grid{2});
        P = zeros(size(Y1));
        for s = 1:M.n_subpop
            switch M.distribution{s,e}
                case {'logn','logn_median','logn_mean'}
                    P = P + w{s}(k)*reshape(bsxfun(@rdivide,mvnpdf(log([Y1(:),Y2(:)]),mu{s}(k,:),permute(Sigma{s}(k,:,:),[2,3,1])),prod([Y1(:),Y2(:)],2)),size(Y1));
                case 'norm'
                    P =  P + w{s}(k)*reshape(mvnpdf([Y1(:),Y2(:)],mu{s}(k,:),permute(Sigma{s}(k,:,:),[2,3,1])),size(Y1));
            end
        end
    end
end
% Plot
if (~options.replicates && ~isempty(D(e).y) && plotData) || ...
        (options.replicates && ~isempty(D(e).replicate(r).y) && plotData)
    if ~isempty(y)
        if D(e).n_dim==1 || ind > 0
            switch options.data.plot
                case 'filled'
                    legendhandles.data = fill(y_hist(round(0.5:0.5:length(y_hist))),[0;h(round(0.5:0.5:length(h)));0],options.data.fill_col{e}); hold on;
                    legendhandles.data.EdgeColor = options.data.fill_col{e};
                case 'empty'
                    legendhandles.data = plot(y_hist(round(0.5:0.5:length(y_hist))),[0;h(round(0.5:0.5:length(h)));0],...
                        '-','color',options.data.col{e},'linewidth',options.data.lw); hold on;
            end
        else
            hs=scatter(log10(y(:,1)),log10(y(:,2)),options.data.markersize,'.'); hold on;
            set(hs,'MarkerEdgeColor',options.data.col{e});
            set(hs,'MarkerEdgeAlpha',0.1);
            [~,kdensity,X1,X2]=kde2d(log10(y));
            contour(X1,X2,kdensity,options.model.levelsets{e,d},'color',options.model.col{e},...
                'LineWidth',options.model.level_linewidth); hold on;
            xlim([log10(y_grid{1}(1)),log10(y_grid{1}(end))]);
            xlim([log10(y_grid{2}(1)),log10(y_grid{2}(end))]);
            box on;
        end
    end
end
if ~isempty(M) && plotModel
    if  D(e).n_dim==1 || ind > 0
        switch options.model.plot
            case 'continuous'
                switch options.x_scale
                    case 'lin'
                        if options.model.subpopulations
                            for s = 1:M.n_subpop
                                plot(y_grid,p_s{s}*d_y_hist,'--','color',...
                                    options.model.col{e},'linewidth',options.model.lw,...
                                    'linestyle',options.model.ls); hold on;
                            end
                        end
                        legendhandles.model = plot(y_grid,p*d_y_hist,'-','color',...
                            options.model.col{e},'linewidth',options.model.lw,  'linestyle',options.model.ls);
                    case 'log'
                        if options.model.subpopulations
                            for s = 1:M.n_subpop
                                legendhandles.subpop = plot(y_grid,p_s{s}.*y_grid*d_y_hist*log(10),'--','color',options.model.col{e},'linewidth',options.model.lw); hold on;
                            end
                        end
                        legendhandles.model = plot(y_grid,p.*y_grid*d_y_hist*log(10),'-',...
                            'color',options.model.col{e},'linewidth',options.model.lw,...
                            'linestyle',options.model.ls);
                end
            case 'hist'
                h_sim = diff(cp)*d_y_hist/d_y_grid;
                plot(y_grid(round(0.5:0.5:length(y_hist))),[0;h_sim(round(0.5:0.5:length(h_sim)));0],...
                    '-','color',options.model.col{e},'linewidth',options.model.lw); hold on;
        end
    else
        if strcmp(options.x_scale,'log')
            P = P.*Y1.*Y2.*log(10)^2;
        end
        hs=scatter(y(:,1),y(:,2),options.data.markersize,'.'); hold on;
        set(hs,'MarkerEdgeColor',options.data.col{e});
        set(hs,'MarkerEdgeAlpha',0.1);
        box on;
        contour(Y1,Y2,P,options.model.levelsets{e,d},'color',options.model.col{e},'LineWidth',options.model.level_linewidth); hold on;
    end
end

if D(e).n_dim == 1 || ind > 0
    
    if ~options.plainstyle
        if ind > 0
            xlabel(D(e).measurand{ind});
        else
            xlabel(D(e).measurand);
        end
        ylabel('frequency');
    end
    if plotModel
        xlim([y_min{e},y_max{e}]);
        if ~isempty(options.z_max{e})
            ylim([0,options.z_max{e}]);
        end
    end
    if legendflag
        if isempty(M)
            legend('data')
        else
            if options.model.subpopulations
                legend([legendhandles.data,legendhandles.model,legendhandles.subpop],...
                    'data','model','subpopulations')
            else
                legend('data','model')
            end
        end
    end
else
    if plotData
        xlim([log10(y_min{e}(1)),log10(y_max{e}(1))]);
        ylim([log10(y_min{e}(2)),log10(y_max{e}(2))]);
    end
    if plotModel
        set(gca,'xscale',options.x_scale);
        set(gca,'yscale',options.x_scale);
    end
    if ~options.sameplot || mod(e,2)
        if plotModel
            xlim([y_min{e}(1),y_max{e}(1)]);
            ylim([y_min{e}(2),y_max{e}(2)])
        end
    else
        xlim([y_min{e-1}(1),y_max{e-1}(1)]);
        ylim([y_min{e-1}(2),y_max{e-1}(2)])
    end
    if options.switch_axes
        set(gca,'ydir','reverse');
    end
    if ~options.plainstyle
        ylabel(D(e).measurand{1});
        xlabel(D(e).measurand{2});
        
        if ~options.data.kde & legendflag
            legend('data','model')
            colormap(options.model.colormap{e})
        end
        
    end
end
if ~isempty(options.xtick{e})
    set(gca,'xtick',options.xtick{e})
end
if ~isempty(options.ytick{e})
    set(gca,'ytick',options.ytick{e})
end
if plotModel || D(e).n_dim == 1
    set(gca,'xscale',options.x_scale);
end

end
