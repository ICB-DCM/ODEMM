clear all


distributions = {'norm','skew_norm','students_t'};

number_of_previously_used_cores = maxNumCompThreads(1);

for iDist = 1:length(distributions)
    clear parameteters D M options conditions
    load('./data/log_data_PDL');
    Dtmp = D;
    clear D
    D(1) = Dtmp(1);
    D(2) = Dtmp(2);
    
    parameters.max = [6*ones(7,1);... %kinetic
        6*ones(3,1);...%scaling
        6*ones(2,1); % noise
        1]; % weight
    parameters.min = [-6*ones(7,1);... %kinetic
        -6*ones(3,1);...%scaling
        -6*ones(2,1); % noise
        0]; %
    
    parameters.name = {'log_{10}(k_1)',...
        'log_{10}(k_2)',...
        'log_{10}(k_4)',...
        'log_{10}(k_5)',...
        'log_{10}(\beta_{k_3[TrkA]_{0,1}})',...
        'log_{10}(\beta_{k_3[TrkA]_{0,2}})',...
        'log_{10}(\beta_{c_P^{1}[Erk]_{0}})',...
        'log_{10}(o_P^{1})',...
        'log_{10}(c_P^{2})',...
        'log_{10}(o_P^{2})',...
        'log_{10}(\sigma_{P,noise}^{1})',...
        'log_{10}(\sigma_{P,noise}^{2})',...
        'w_1'};
    
    if strcmp(distributions{iDist},'skew_norm')
        parameters.name{14} = '\\delta_P';
        parameters.max(14) = 2;
        parameters.min(14) = -2;
    end
    
    M.n_subpop = 2; %number of subpopulations
    M.model = @(T,theta,u)simulate_SigmaPoints_sPsET_loglog_corr(T,theta,u(1));
    for s = 1:M.n_subpop
        for e=1:2
            M.mean_ind{s,e} = [1]; %index of mean in output
            M.var_ind{s,e} = [4]; %index of variances in output (empty if RREs used)
        end
        
    end
    
    M.sim_type = 'HO';
    options.replicates = false;
    options.write_parameters = true;
    options.measurement_noise = true;
    options.noise_model = 'additive';
    parameters.number = length(parameters.name);
    
    % definition of symbolic parameters required for generateODEMM
    xi = sym(zeros(parameters.number,1));
    for i = 1:parameters.number
        xi(i) = sym(['xi_' num2str(i,'%d')]);
    end
    u = sym(zeros(2,1)); % 1 for definition of stimulation, 1 for subpop,
    for i = 1:7
        u(i) = sym(['u_' num2str(i,'%d')]);
    end
    g = size(D(e).u,1); % number of input for simulation of inidividual subpopulation
    
    for e=[1:2]
        M.u{1,e} = [0]; % subpop
        M.u{2,e} = [1];
    end
    M.sym.theta = [log(10)*(xi(1));... %k1
        log(10)*(xi(2));... %k2
        log(10)*(xi(3));... %k4
        log(10)*(xi(4));... %k5
        log(10)*(u(g+1)*xi(5) + (1-u(g+1))*xi(6));... %k3TrkA
        log(10)*(xi(7));...  %sPmErk
        log(log(10.^(2*(u(g+1)*xi(8) + (1-u(g+1))*xi(9)))+1));... % log(sigma^2)
        log(10)*xi(11);...
        log(log(10.^(2*(xi(10)))+1))];
    
    r=1;
    e=1;
    M.sym.scaling{r,e} = sym('1');
    M.sym.offset{r,e} = 10.^(xi(8));
    M.sym.sigma_noise{e} = 10.^(xi(11));
    if strcmp(distributions{iDist},'skew_norm')
        M.sym.delta{1,e} = xi(14);
        M.sym.delta{2,e} = xi(14);
    end
    
    e=2;
    M.sym.scaling{r,e} = 10.^(xi(9));
    M.sym.offset{r,e} = 10.^(xi(10));
    M.sym.sigma_noise{e} = 10.^(xi(12));
    if strcmp(distributions{iDist},'skew_norm')
        M.sym.delta{1,e} = xi(14);
        M.sym.delta{2,e} = xi(14);
    end
    
    for e = 1:2
        for s = 1:2
            M.distribution{s,e} = distributions{iDist};
        end
        M.sym.w{1,e} = ones(size(D(e).t(:)))*xi(13);
        M.sym.w{2,e} = ones(size(D(e).t(:)))*(1-xi(13));
    end
    M.name = ['NGFErk_1D_' distributions{iDist}];
    options.simulate_musigma = false;
    [conditions,D] = collectConditions(D,M);
    generateODEMM_extended(D,M,parameters,conditions,options);
    eval(['ODEMM_' M.name])
    
    options.MS = PestoOptions();
    options.MS.n_starts = 100;
    options.MS.localOptimizerOptions.Display = 'iter';
    options.MS.foldername = ['./results/results_' M.name];
    options.MS.save = true;
    warning off
    parameters.guess = getParameterGuesses(parameters,@(xi) ...
        logLikelihood_extend(xi,M,D,options,conditions),...
        options.MS.n_starts, parameters.min,parameters.max);
    warning off
    
    parameters = getMultiStarts(parameters,@(xi) ...
        logLikelihood_extend(xi,M,D,options,conditions),options.MS);
    save(options.MS.foldername)
end