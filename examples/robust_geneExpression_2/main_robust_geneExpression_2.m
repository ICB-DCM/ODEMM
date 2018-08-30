clear all
close all
clc

load data/data_geneExpression_2_SSA
n_data = sum(sum(~isnan(D(1).y)));

distributions = {'norm','skew_norm','students_t','neg_binomial'};%,'logn_mean'};

for iDist = [1,3,4]%1:length(distributions)
    clear parameters M
    if strcmp(distributions{iDist},'skew_norm')
        parameters.name = {'log_{10}(k_{1})',...
            'log_{10}(k_{2,1})',...
            'log_{10}(k_{2,2})',...
            'log_{10}(k_3)',...
            'log_{10}(k_4)',...
            'w',...
            '\\delta'};
    else
        parameters.name = {'log_{10}(k_{1})',...
            'log_{10}(k_{2,1})',...
            'log_{10}(k_{2,2})',...
            'log_{10}(k_3)',...
            'log_{10}(k_4)',...
            'w'};
    end
    parameters.number = length(parameters.name);
    M.n_subpop = 2; %number of subpopulations
    M.model = @(T,theta,u) simulate_geneExpression_2_MA(T,theta,[]);
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
    M.sym.theta = [10.^xi(1);u(n_siminput+1)*10.^xi(2)+(1-u(n_siminput+1))*10.^xi(3);10.^xi(4);10.^xi(5)];
    
    if strcmp(distributions{iDist},'skew_norm')
        M.sym.delta{1,e} = xi(7);
        M.sym.delta{2,e} = xi(7);
    end
    M.u{1,e} = [1]; % subpop 1
    M.u{2,e} = [0]; % subpop 2
    % initialize weights
    M.sym.w{1,e} = xi(6)*ones(size(D(e).t(:)));
    M.sym.w{2,e} = (1-xi(6))*ones(size(D(e).t(:)));
    r=1; % only one replicate
    M.sym.scaling{r,e} = sym('1');
    M.sym.offset{r,e} = sym('0');
    
    
    %% distribution assumption
    for s = 1:M.n_subpop
        M.distribution{s,e} = distributions{iDist};
    end
    
    [conditions,D] = collectConditions(D,M);
    options.dimension = 'univariate'; % 1D-measurements
    
    parameters.min = [-3*ones(5,1);0];
    parameters.max = [ 3*ones(5,1);1];
    if strcmp(distributions{iDist},'skew_norm')
        parameters.min(end+1) = -5e2;
        parameters.max(end+1) = 5e2;
    end
    
    M.name = [M.distribution{1,1} '_geneExpression_2'];
    
    generateODEMM_extended(D,M,parameters,conditions,options);
    eval(['ODEMM_' M.name]);
    
    if strcmp(M.distribution{s,e},'students_t')
        parameters.max(end) = 8;
    end
    %     xi = rand(parameters.number,1);
    %     ll = logLikelihood_extend(xi,M,D,options);
    %     [ll,dll] =  logLikelihood_extend(xi,M,D,options);
    %     [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logLikelihood_extend(xi,M,D,options),1e-6);
    %     [g,g_fd_f,g_fd_b,g_fd_c]
    
    options.MS = PestoOptions();
    options.MS.n_starts = 30;
    options.MS.localOptimizerOptions.Display = 'iter';
    warning off

    parameters.guess = getParameterGuesses(parameters,@(xi) ...
        logLikelihood_extend(xi,M,D,options,conditions),...
        options.MS.n_starts, parameters.min,parameters.max);
    warning off
%     if strcmp(distributions{iDist},'skew_norm')
%         @(xi) confuneq_skew_norm(xi,M,D,conditions)
%     end
    load(['./results/results_geneExpression_2_SSA_' distributions{iDist}])

    %parameters = getMultiStarts(parameters,@(xi) ...
       %logLikelihood_extend(xi,M,D,options,conditions),options.MS);
    parameters = getParameterProfiles(parameters,@(xi) ...
       logLikelihood_extend(xi,M,D,options,conditions),options.MS);
    xi = parameters.MS.par(:,1);
    plotODEMM(D,M,xi)
    save(['./results/results_geneExpression_2_SSA_' distributions{iDist}])
    
%     for jDist = 1:length(distributions)
%         load(['data/data_geneExpression_2_MA_' distributions{jDist}])
%         [conditions,D] = collectConditions(D,M);
%         parameters = getMultiStarts(parameters,@(xi) ...
%             logLikelihood_extend(xi,M,D,options,conditions),options.MS);
%         save(['./results/results_geneExpression_2_MA_' ...
%                 distributions{jDist} '_' distributions{iDist}])
%     end
end
