function lppd = run_lppd_cr(model_name,burnin,nIter,saveFlag)
% This function runs the sampling for the log pointwise predictive density
% for the conversion reaction reaction example
%
% USAGE:
% parameters_lppd = run_lppd_cr(model_name,burnin,nIter,saveFlag)
%
% Parameters:
%   model_name: string indicating model ('RRE_timedep','RRE_subpop',
%               'RRE_onlyone','SP_k1','SP_k2','SP_k3','SP_k1k2','SP_k1k3',
%               'SP_k2k3',or 'SP_all')
%   burnin: number of burn-in steps for sampling
%   nIter: number of total iterations for sampling
%   saveFlag: if results should be saved
%
% Return values:
%   lppd: log pointwise predictive density

load(['./results/results_' model_name])
if isfield(options,'measurement_noise')
    options.llh.measurement_noise = options.measurement_noise;
    options.llh.noise_model = options.noise_model;
    options.llh.simulate_musigma = options.simulate_musigma;
else
    options.llh.dimension = options.dimension;
    options.llh.sigmas = options.sigmas;
end


options.llh.prior.flag = true;
options.llh.prior.distribution = 'uniform';
options.llh.prior.min = parameters.min;
options.llh.prior.max = parameters.max;
options.pesto = PestoOptions();

% Split data set in 5 sets
for iData = 1:5
    D_train{iData} = D;
    D_train{iData}.y(1,iData,:) = nan;
    
    D_test{iData} = D;
    for i = 1:5
        if ~isequal(iData,i)
            D_test{iData}.y(1,i,:) = nan;
        end
    end
end

% Sample for each splitted data set
for iData = 1:5
    options.pesto.MCMC.nIterations         = nIter;
    options.pesto.MCMC.sigma0              = 10*diag(ones(1,parameters.number));
    options.pesto.MCMC.theta0              = parameters.MS.par(:,1); % starting in the MLE
    options.pesto.MCMC.samplingAlgorithm   = 'PT';
    options.pesto.MCMC.PT.nTemps           = 10;
    options.pesto.mode                     = 'text';
    timestart = tic;
    parameters_lppd{iData} = getParameterSamples(parameters, ...
        @(xi) logLikelihood(xi,M,D_train{iData},options.llh,conditions), options.pesto);
    
    [z,~]=gewekeTest(squeeze(parameters_lppd{iData}.S.par(:,burnin+1:end,1))',0.1,0.5);
    parameters_lppd{iData}.S.max_zscore = max(abs(z));
    parameters_lppd{iData}.S.burnin = burnin;
    parameters_lppd{iData}.S.t_cpu = toc(timestart);
    if saveFlag
        save(['results_lppd_' model_name])
    end
end

% Calculate log pointwise predictive density
options.llh.prior.flag = false;
lppd = 0;

% for each set of the data evaluate prediction
for iData = 1:5
    clear q_i
    iCount = 1;
    
    % thinning such that 5000 data points are used
    thin = ceil(numel(parameters_lppd{iData}.S.burnin+1:size(parameters_lppd{iData}.S.par,2))/5000);
    for iSample = parameters_lppd{iData}.S.burnin+1:thin:size(parameters_lppd{iData}.S.par,2)
        q_i(iCount) = logLikelihood(parameters_lppd{iData}.S.par(:,iSample,1),...
            M,D_test{iData},options.llh,conditions);
        iCount = iCount+1;
    end
    w = 1./sum(~isnan(q_i))*ones(size(q_i));
    
    lppd = lppd + computeMixtureProbability(w,q_i);
end

% save results to file
if saveFlag
    save(['./results/results_lppd_' model_name],'parameters_lppd','lppd',...
        'options','M','D','conditions','model_name')
end
end


