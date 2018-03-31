function Q = run_thermIntegration_cr(model_name,burnin,nIter,saveFlag)
% This function runs thermodynamic integration for the conversion reaction
% example.
%
% USAGE:
% Q = run_thermIntegration_cr(model_name,burnin,nIter,saveFlag)
%
% Parameters:
%   model_name: string indicating model ('RRE_timedep','RRE_subpop','RRE_onlyone','SP_k1','SP_k2',
%        'SP_k3','SP_k1k2','SP_k1k3','SP_k2k3',or 'SP_all')
%   burnin: number of burn-in steps for sampling
%   nIter: number of total iterations for sampling
%   saveFlag: if results should be saved
%
% Return values:
%   Q: log marginal likelihood

load(['./results/results_' model_name])
if isfield(options,'measurement_noise')
    options.llh.measurement_noise = options.measurement_noise;
    options.llh.noise_model = options.noise_model;
    options.llh.simulate_musigma = options.simulate_musigma;
else
    options.llh.dimension = options.dimension;
    options.llh.sigmas = options.sigmas;
end

%% initialize log posterior
options.llh.prior.flag = true;
options.llh.prior.distribution = 'uniform';
options.llh.prior.min = parameters.min;
options.llh.prior.max = parameters.max;

q = 5;
phi=@(x) x.^q;
phiprime = @(x) q*(x.^(q-1));

%% perform sampling for different temperatures
options.llh.prior.flag = true; % compute only likelihoood
nTemps = 17;
lambda = linspace(0,1,nTemps); % these are equidistant
temperatures = phi(lambda); %these are not equidistant
timestart1 = tic;
options.pesto = PestoOptions();

for l = 1:numel(temperatures)
    options.pesto.MCMC.nIterations         = nIter;
    options.pesto.MCMC.sigma0              = 10*diag(ones(1,parameters.number));
    options.pesto.MCMC.theta0 = parameters.MS.par(:,1);
    options.pesto.MCMC.samplingAlgorithm   = 'PT';
    options.pesto.MCMC.PT.nTemps           = 7;
    options.pesto.mode = 'text';
    
    options.llh.tau = temperatures(l);
    timestart = tic;
    
    parameters_l{l} = getParameterSamples(parameters, ...
        @(xi) logLikelihood(xi,M,D,options.llh,conditions), options.pesto);
    [z,~]=gewekeTest(squeeze(parameters_l{l}.S.par(:,burnin+1:end,1))',0.1,0.5);
    parameters_l{l}.S.max_zscore = max(abs(z));
    parameters_l{l}.S.burnin = burnin;
    parameters_l{l}.S.t_cpu = toc(timestart);

    if saveFlag
        save(['./results/results_BFchains_' model_name])
    end
end
if saveFlag
    save(['./results/results_BFchains_' model_name])
end

% calculate expected log deviance for every temperature
for iTemp = 1:numel(temperatures)
    sum_ll = 0;
    thin = ceil(numel(parameters_l{iTemp}.S.burnin+1:...
        size(parameters_l{iTemp}.S.par,2))/5000);
    for j=parameters_l{iTemp}.S.burnin+1:thin:size(parameters_l{iTemp}.S.par,2)
        ll = logLikelihood(parameters_l{iTemp}.S.par(:,j),M,D,options.llh,conditions);
        sum_ll= sum_ll+ll;
    end
    f(iTemp) = phiprime(lambda(iTemp))*sum_ll/...
        numel([parameters_l{iTemp}.S.burnin+1:thin:...
        size(parameters_l{iTemp}.S.par,2)]);
end
clear iTemp ll sum_ll j thin 

% Simpsons rule to calculate integral
h = lambda(3) - lambda(1);
Q = (h/6)*(f(1) + f(end) + sum(2*f(3:2:end-2)) + 4*sum(f(2:2:end-1)));
timelogdev = toc(timestart1);

if saveFlag
    save(['./results/results_BFchains_' model_name])
end
end
