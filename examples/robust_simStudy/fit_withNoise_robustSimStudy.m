clear all
close all
clc

outlierflag = false;
outlierstr = 'outlier5_dublets'; % 'outlier10_unif', 'outlier2_zeros' 
% (only used if outlierflag == true)

load_simStudy_settings

modelnames = {'conversionReaction','diffProteinExpression',...
    'twoStageGeneExpression','diffProteinExpression'};
distributions = {'norm','skew_norm','students_t','neg_binomial'};

for m = [1,3,4]
    for it = 1:3
        t = tps{it};
        for ic = 1:4
            for set = 1:3
                if outlierflag
                    load(['./dataOutlier/data_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                        'cells_' num2str(length(t)) 'tps_' num2str(set) 'paramset_' outlierstr],'D')
                else
                    load(['./data/data_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                        'cells_' num2str(length(t)) 'tps_' num2str(set) 'paramsetD'],'D')
                end
                D(1).u = 1;
                nDist = 4;
                for iDist = 4
                    if outlierflag
                        checkExist = exist(['./resultsOutlierNoise/results_noise_model' ...
                            num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                            'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                            'paramset_' distributions{iDist} '_' outlierstr '.mat'],'file');
                    else
                        checkExist = exist(['./resultsNoise/results_noise_model' ...
                            num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                            'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                            'paramset_' distributions{iDist} '.mat'],'file');
                    end
                    if  checkExist == 0
                        try
                            clear parameters M
                            switch m
                                case 1
                                    parameters.name = {'log_{10}(k_{1,1})',...
                                        'log_{10}(k_{1,2})',...
                                        'log_{10}(k_{2})',...
                                        'log_{10}(k_3)',...
                                        'w',...
                                        '\sigma_{noise}'};
                                    ind_weight = 5;
                                    ind_noise = 6;
                                case 2
                                    error('adapt')
                                    parameters.name = {'log_{10}(lambda0)',...
                                        'log_{10}(lambdaA_{1})',...
                                        'log_{10}(lambdaA_{2})',...
                                        'log_{10}(lambdaB_{1})',...
                                        'log_{10}(lambdaB_{2})',...
                                        'log_{10}(gammaAB)',...
                                        'w',...
                                        '\sigma_{noise}'};
                                    ind_weight = 7;
                                    ind_noise = 8;
                                case 3
                                    parameters.name = {'log_{10}(k_{1})',...
                                        'log_{10}(k_{2,1})',...
                                        'log_{10}(k_{2,2})',...
                                        'log_{10}(k_3)',...
                                        'log_{10}(k_4)',...
                                        'log_{10}(k_5)',...
                                        'w',...
                                        '\sigma_{noise}'};
                                    ind_weight = 7;
                                    ind_noise = 8;
                                case 4
                                    parameters.name = {'log_{10}(lambda0)',...
                                        'log_{10}(lambdaA_{1})',...
                                        'log_{10}(lambdaA_{2})',...
                                        'log_{10}(gammaAB)',...
                                        'w',...
                                        '\sigma_{noise}'};
                                    ind_weight = 5;
                                    ind_noise = 6;
                            end
                            
                            if strcmp(distributions{iDist},'skew_norm')
                                if m == 2
                                    parameters.name{end+1} = '\\delta_A';
                                    parameters.name{end+1} = '\\delta_B';
                                else
                                    parameters.name{end+1} = '\\delta';
                                end
                            end
                            parameters.number = length(parameters.name);
                            M.n_subpop = 2; %number of subpopulations
                            eval(['M.model = @(T,theta,u) simulate_model' num2str(m) '_' modelnames{m} '(T,theta,[])']);
                            
                            e=1; % experiment index (only one experiment)
                            if m == 1 || m == 3 || m == 4
                                for s = 1:2
                                    M.mean_ind{s,e} = [1]; %index of mean in output
                                    M.var_ind{s,e} = [2]; %index of variances in output (empty if RREs used)
                                end
                            else
                                for s = 1:2
                                    M.mean_ind{s,e} = [1,2]; %index of mean in output
                                    M.var_ind{s,e} = [3,4,5]; %index of variances in output (empty if RREs used)
                                end
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
                            switch m
                                case 1
                                    M.sym.theta = [10.^(u(2)*xi(1)+(1-u(2))*xi(2));...
                                        10.^xi(3);...
                                        10.^xi(4)];
                                case 2
                                    M.sym.theta = [10.^xi(1);...
                                        10.^(u(2)*xi(2)+(1-u(2))*xi(3));...
                                        10.^(u(2)*xi(4)+(1-u(2))*xi(5));
                                        10.^xi(6)];
                                case 3
                                    M.sym.theta = [10.^xi(1);
                                        10.^(u(2)*xi(2)+(1-u(2))*xi(3));...
                                        10.^xi(4);...
                                        10.^xi(5);
                                        10.^xi(6)];
                                case 4
                                    M.sym.theta = [10.^xi(1);...
                                        10.^(u(2)*xi(2)+(1-u(2))*xi(3));...
                                        10.^xi(4)];
                            end
                            if strcmp(distributions{iDist},'skew_norm')
                                if m == 2
                                    M.sym.delta{1,e} = xi(parameters.number-1:end);
                                    M.sym.delta{2,e} = xi(parameters.number-1:end);
                                else
                                    M.sym.delta{1,e} = xi(parameters.number);
                                    M.sym.delta{2,e} = xi(parameters.number);
                                end
                            end
                            
                            M.u{1,e} = [1]; % subpop 1
                            M.u{2,e} = [0]; % subpop 2
                            
                            % initialize weights
                            M.sym.w{1,e} = xi(ind_weight)*ones(size(D(e).t(:)));
                            M.sym.w{2,e} = (1-xi(ind_weight))*ones(size(D(e).t(:)));
                            r=1; % only one replicate
                            if m == 2
                                M.sym.scaling{r,e} = [sym('1');sym('1')];
                                M.sym.offset{r,e} = [sym('0');sym('0')];
                            else
                                M.sym.scaling{r,e} = sym('1');
                                M.sym.offset{r,e} = sym('0');
                            end
                            M.sym.sigma_noise{e} = [10.^xi(ind_noise)];
                            
                            
                            %% distribution assumption
                            for s = 1:M.n_subpop
                                M.distribution{s,e} = distributions{iDist};
                            end
                            
                            [conditions,D] = collectConditions(D,M);
                            options.dimension = 'univariate'; % 1D-measurements
                            
                            parameters.min = -3*ones(parameters.number,1);
                            parameters.min(ind_weight) = 0;
                            parameters.max = 3*ones(parameters.number,1);
                            parameters.max(ind_weight) = 1;
                            
                            if strcmp(distributions{iDist},'skew_norm')
                                parameters.min(parameters.number) = -5e2;
                                parameters.max(parameters.number) = 5e2;
                            end
                            if outlierflag
                                M.name = ['noise_' M.distribution{1,1} '_' modelnames{m} '_' outlierstr];
                                if m == 3 && iDist == 4 && strcmp(outlierstr,'outlier5_dublets')
                                    M.name = ['noise_' M.distribution{1,1} '_twoStage_' outlierstr];
                                end
                            else
                                M.name = ['noise_' M.distribution{1,1} '_' modelnames{m}];
                            end
                            options.measurement_noise = true;
                            options.noise_model = 'additive';
                            generateODEMM_extended(D,M,parameters,conditions,options);
                            eval(['ODEMM_' M.name ';']);
                            
                            if strcmp(M.distribution{s,e},'students_t')
                                parameters.max(end) = 8;
                            end
                            %xi = rand(parameters.number,1);
                            %     ll = logLikelihood_extend(xi,M,D,options);
                            %     [ll,dll] =  logLikelihood_extend(xi,M,D,options);
                            %[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logLikelihood_extend(xi,M,D,options),1e-6);
                            %[g,g_fd_f,g_fd_b,g_fd_c]
                            
                            options.MS = PestoOptions();
                            options.MS.n_starts = 30;
                            options.MS.localOptimizerOptions.Display = 'iter';
                            warning off
                            
                            parameters.guess = getParameterGuesses(parameters,@(xi) ...
                                logLikelihood_extend(xi,M,D,options,conditions),...
                                options.MS.n_starts, parameters.min,parameters.max);
                            warning off
                            parameters = getMultiStarts(parameters,@(xi) ...
                                logLikelihood_extend(xi,M,D,options,conditions),options.MS);
                            %parameters = getParameterProfiles(parameters,@(xi) ...
                            %    logLikelihood_extend(xi,M,D,options,conditions),options.MS);
                            if outlierflag
                                save(['./resultsOutlierNoise/results_noise_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                                    'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                                    'paramset_' distributions{iDist} '_' outlierstr]);
                            else
                                save(['./resultsNoise/results_noise_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                                    'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                                    'paramset_' distributions{iDist}]);
                            end
                        catch e
                            fprintf(1,'error  message:\n%s',e.message);
                        end
                    else
                        disp(['model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                            'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                            'paramset_' distributions{iDist} ' skipped']);
                    end
                end
            end
        end
    end
end