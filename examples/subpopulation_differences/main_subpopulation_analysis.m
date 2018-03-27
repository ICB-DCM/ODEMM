clear all
close all
clc

% Generate the model files for all combinations of differences
generate_nosubpop_file % is icomb = 1 with no difference
generate_subpop_files % icomb 2,...,64


%% Optimization
% Run optimization for each model, note that it will take quite long and
% was parallelized for the results in the paper

for icomb = 1:64
    parameters{icomb} = run_fitting_subpop_icomb(icomb);
end

save ./results/results_subpop_differences parameters

%% Optimize final model with different coefficient of variations 
generate_subpop_TrkA

load('./data/data_PDL');
options.MS = PestoOptions();
options.MS.n_starts = 100;
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('GradObj','on',...
    'display','iter','TolFun',1e-10, 'TolX',1e-14, 'MaxIter', 1000,...
    'algorithm','interior-point');
options.MS.comp_type = 'sequential';
options.MS.mode = 'text';
options.simulate_musigma = true;

ODEMM_NGF_subpop_TrkA

options.MS.foldername = ['results_' M.name];
[conditions,D] = collectConditions(D,M);

% Multi-start optimiization
parameters.guess = getParameterGuesses(parameters,@(xi) ...
    logLikelihood(xi,M,D,options,conditions),...
    options.MS.n_starts,parameters.min,parameters.max);
parameters = getMultiStarts(parameters,@(xi) ...
    logLikelihood(xi,M,D,options,conditions),options.MS);

% Profile likelihood
parameters = getParameterProfiles(parameters,@(xi) ...
    logLikelihood(xi,M,D,options,conditions),options.MS);

%% Sampling: lppd 
count = 1;
for icomb = [1,2,3,6,9,17,33]
    llpds(count) = run_lppd_subpop(icomb);
end

save ./results/results_subpop_lppds lppds

%% Sampling: log-marginals/Bayes factor
count = 1;
for icomb = [1,2,3,6,9,17,33]
    logmargs(count) = run_logmarg_subpop(icomb);
end

save ./results/results_subpop_logmargs logmargs
