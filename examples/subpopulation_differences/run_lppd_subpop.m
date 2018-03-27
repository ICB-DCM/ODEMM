function lppd = run_lppd_subpop(icomb)
% This function calculates the log pointwise predictive density for 
% the model of NGF-induced Erk1/2 signaling accounting for the differences 
% in subpopulations denoted by icomb. 

load('./data/data_PDL');
load('./results/results_subpop_differences');
ptemp = parameters{icomb};
clear parameters
parameters = ptemp;

eval(['ODEMM_NGF_subpop_comb' num2str(icomb)]);
[conditions,D] = collectConditions(D,M);

options.llh.simulate_musigma = true;
options.llh.prior.flag = true;
options.llh.prior.distribution = 'uniform';
options.llh.prior.min = parameters.min;
options.llh.prior.max = parameters.max;

iData = 1;
D_train{iData} = D;
D_test{iData} = D;
e=1; %kinetic
D_train{iData}(e).y(1,[2,3,5,6,7],:) = nan;
D_test{iData}(e).y(1,[1,4],:) = nan;
e=2; %kinetic
D_train{iData}(e).y([2,3,5,6,7],1,:) = nan;
D_test{iData}(e).y([1,4],1,:) = nan;
for e=3:4
    D_train{iData}(e).y([2,3,5,6],1,:,:) = nan;
    D_test{iData}(e).y([1,4],1,:,:) = nan;
end

iData = 2;
D_train{iData} = D;
D_test{iData} = D;
e=1; % kinetic
D_train{iData}(e).y(1,[1,3,4,6,7],:) = nan;
D_test{iData}(e).y(1,[2,5],:) = nan;
e=2; % dose response
D_train{iData}(e).y([1,3,4,6,7],1,:) = nan;
D_test{iData}(e).y([2,5],1,:) = nan;
for e=3:4 % dose response
    D_train{iData}(e).y([1,3,4,6],1,:,:) = nan;
    D_test{iData}(e).y([2,5],1,:,:) = nan;
end

iData = 3;
D_train{iData} = D;
D_test{iData} = D;
e=1; % kinetic
D_train{iData}(e).y(1,[1,2,4,5],:) = nan;
D_test{iData}(e).y(1,[3,6,7],:) = nan;
e=2; % dose response
D_train{iData}(e).y([1,2,4,5,],1,:) = nan;
D_test{iData}(e).y([3,6,7],1,:) = nan;
for e=3:4 % dose response
    D_train{iData}(e).y([1,2,4,5],1,:,:) = nan;
    D_test{iData}(e).y([3,6],1,:,:) = nan;
end

burnin = 1e5;
nIter = 3e5;
options.pesto = PestoOptions();
for iData = 1:3
    options.pesto.MCMC.nIterations         = nIter;
    options.pesto.MCMC.sigma0              = 1e-5*diag(ones(1,parameters.number));
    options.pesto.MCMC.theta0 = parameters.MS.par(:,1);
    options.pesto.MCMC.samplingAlgorithm   = 'PT';
    options.pesto.MCMC.PT.nTemps           = 7;
    options.pesto.mode = 'text';
    timestart = tic;
    parameters_lppd{iData} = getParameterSamples(parameters, ...
        @(xi) logLikelihood(xi,M,D_train{iData},options.llh,conditions), options.pesto);
    [z,~]=gewekeTest(squeeze(parameters_lppd{iData}.S.par(:,burnin+1:end,1))',0.1,0.5);
    parameters_lppd{iData}.S.max_zscore = max(abs(z));
    parameters_lppd{iData}.S.burnin = burnin;
    parameters_lppd{iData}.S.t_cpu = toc(timestart);
end

save(['results_subpop_lppd_' num2str(icomb)])
   
options.llh.prior.min = nan;
options.llh.prior.max = nan;
lppd = 0;
for iData = 1:3
    iCount = 1;
    options.llh.prior.flag = 0;
    thin = ceil(numel(parameters_lppd{iData}.S.burnin+1:size(parameters_lppd{iData}.S.par,2))/1000);
    for iSample = parameters_lppd{iData}.S.burnin+1:thin:size(parameters_lppd{iData}.S.par,2)
        q_i(iCount) = logLikelihood(parameters_lppd{iData}.S.par(:,iSample,1),...
            M,D_test{iData},options.llh,conditions);
        iCount = iCount+1;
    end
    w = 1./numel(q_i)*ones(size(q_i));
    lppd = lppd + computeMixtureProbability(w,q_i);
end

save(['./results/results_subpop_lppd_' num2str(icomb)])
    
end
