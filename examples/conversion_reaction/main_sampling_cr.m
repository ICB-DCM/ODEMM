% This is the main script for the Bayesian model selection for the
% conversion reaction.

clear all
close all
clc

models = {'RRE_onlyone','RRE_subpop','RRE_timedep','SP_k1', 'SP_k1k2',...
    'SP_all', 'SP_k2', 'SP_k2k3', 'SP_k1k3', 'SP_k3'};

saveFlag = true;
nIter = 2e4;
burnin = 1e3; % BUT check results whether chain really converged, Geweke test 
% is often not conservative enough

num_repeats = 1; % how often Bayes factor and lppd should be evaluated for each model

% initialization
logmargs = nan(num_repeats,10);
lppds = nan(num_repeats,10);

for i = 1:num_repeats %repeat sampling
    for m = 1:numel(model_names)
        % log marginal likelihood using thermodynamic integration
        logmargs(i,m) = run_thermIntegration_cr(models{m},burnin,nIter,saveFlag);

        % log pointwise predictive density
        lppds(i,m) = run_lppd_cr(models{m},burnin,nIter,saveFlag);
    end
end

save ./results/results_lppds lppds
save ./results/results_logmargs logmargs

