clear all
close all
clc

load data/data_conversionProcess_SSA
n_data = sum(sum(~isnan(D(1).y)));
distributions = {'students_t','skew_norm','norm','neg_binomial'};

for iDist = 1:length(distributions)
    clear parameters M
    if strcmp(distributions{iDist},'skew_norm')
        parameters.name = {'log_{10}(k_{1,1})',...
            'log_{10}(k_{1,2})',...
            'log_{10}(k_{2})',...
            'log_{10}(k_3)',...
            'w',...
            '\\delta'};
    else
        parameters.name = {'log_{10}(k_{1,1})',...
            'log_{10}(k_{1,2})',...
            'log_{10}(k_{2})',...
            'log_{10}(k_3)',...
            'w'};
    end
    parameters.number = length(parameters.name);
    M.n_subpop = 2; %number of subpopulations
    M.model = @(T,theta,u) simulate_conversionProcess_MA(T,theta,[]);
    e=1; % experiment index (only one experiment)
    for s = 1:2
        M.mean_ind{s,e} = [1]; %index of mean in output
        M.var_ind{s,e} = [2]; %index of variances in output (empty if RREs used)
    end
    M.sim_type = 'HO';
    xi = sym(zeros(parameters.number,1));
    for i = 1:parameters.number
        xi(i) = sym(['xi_' num2str(i,'%d')]);
    end
    n_siminput = size(D(1).u,1); % number of input for simulation of inidividual subpopulation
    u = sym(zeros(n_siminput+1,1));
    for i = 1:n_siminput+1
        u(i) = sym(['u_' num2str(i,'%d')]);
    end
    M.sym.theta = [u(n_siminput+1)*10.^xi(1)+(1-u(n_siminput+1))*10.^xi(2);...
        10.^xi(3);...
        10.^xi(4);];
    if strcmp(distributions{iDist},'skew_norm')
        M.sym.delta{1,e} = 10.^xi(6)-1e3;
        M.sym.delta{2,e} = 10.^xi(6)-1e3;
    end
    M.u{1,e} = [1]; % subpop 1
    M.u{2,e} = [0]; % subpop 2
    % initialize weights
    M.sym.w{1,e} = xi(5)*ones(size(D(e).t(:)));
    M.sym.w{2,e} = (1-xi(5))*ones(size(D(e).t(:)));
    r=1; % only one replicate
    M.sym.scaling{r,e} = sym('1');
    M.sym.offset{r,e} = sym('0');
    
    
    %% distribution assumption
    for s = 1:M.n_subpop
        M.distribution{s,e} = distributions{iDist};
    end
    
    [conditions,D] = collectConditions(D,M);
    options.dimension = 'univariate'; % 1D-measurements
    
    parameters.min = [-3*ones(4,1);0];
    parameters.max = [ 3*ones(4,1);1];
    if strcmp(distributions{iDist},'skew_norm')
        parameters.min(end+1) = 0;
        parameters.max(end+1) = 3.3;
    end
    
    M.name = [M.distribution{1,1} '_conversionProcess'];
    
    generateODEMM_extended(D,M,parameters,conditions,options);
    eval(['ODEMM_' M.name]);
    
    if strcmp(M.distribution{s,e},'students_t')
        parameters.max(end) = 8;
    end
   
    xi = rand(parameters.number,1);
    ll = logLikelihood_extend(xi,M,D,options);
    [ll,dll] =  logLikelihood_extend(xi,M,D,options);
    [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logLikelihood_extend(xi,M,D,options),1e-6);
    [g,g_fd_f,g_fd_b,g_fd_c]
    
    options.MS = PestoOptions();
    options.MS.n_starts = 30;
    warning off
    parameters.guess =   getParameterGuesses(parameters,@(xi) logLikelihood_extend(xi,M,D,options,conditions),...
        options.MS.n_starts, parameters.min,parameters.max);
    warning off
    parameters = getMultiStarts(parameters,@(xi) logLikelihood_extend([xi],M,D,options,conditions),options.MS);
    xi = parameters.MS.par(:,1);
    plotODEMM(D,M,xi)
    save(['./results/results_conversionProcess_SSA_' distributions{iDist}])
    params{iDist}=parameters;
end

for iDist = 1:length(distributions)
    BICs(iDist) = -2*params{iDist}.MS.logPost(1) + log(n_data)*params{iDist}.number;
    llhs(iDist) = params{iDist}.MS.logPost(1);
end

%% Visualization
figure
subplot(2,1,1);
plot(1:length(distributions), BICs-min(BICs), 'o');
set(gca,'xtick',[1:length(distributions)],'xticklabels',distributions)
ylabel('\Delta BIC')
subplot(2,1,2);
plot(1:length(distributions), llhs, 'o')
set(gca,'xtick',[1:length(distributions)],'xticklabels',distributions)
ylabel('log-likelihood')

        