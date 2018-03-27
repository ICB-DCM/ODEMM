function logmarg = run_logmarg_subpop(icomb)
% This function calculates the log marginals for the model of NGF-induced
% Erk1/2 signaling accounting for the differences in subpopulations denoted
% by icomb. 

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
nTemps = 17;
burnin = 1e5;
nIter = 3e5;
q = 5;
phi=@(x) x.^q;
phiprime = @(x) q*(x.^(q-1));
lambda = linspace(0,1,nTemps); % these are equidistant
temperatures = phi(lambda); %these are not equidistant
options.pesto = PestoOptions();
for l = 1:nTemps
    options.pesto.MCMC.nIterations         = nIter;
    options.pesto.MCMC.sigma0              = 1e-3*diag(ones(1,parameters.number));
    options.pesto.MCMC.theta0 = parameters.MS.par(:,1);
    options.pesto.MCMC.samplingAlgorithm   = 'PT';
    options.pesto.MCMC.PT.nTemps           = 5;
    rng(1);
    options.llh.tau = temperatures(l);
    timestart = tic;
    parameters_l{l} = getParameterSamples(parameters, ...
        @(xi) logLikelihood(xi,M,D,options.llh,conditions), options.pesto);
    [z,~]=gewekeTest(squeeze(parameters_l{l}.S.par(:,burnin+1:end,1))',0.1,0.5);
    parameters_l{l}.S.max_zscore = max(abs(z));    
    parameters_l{l}.S.burnin = burnin;
    parameters_l{l}.S.t_cpu = toc(timestart);
end
    
q = 5;
phi=@(x) x.^q;
phiprime = @(x) q*(x.^(q-1));
lambda = linspace(0,1,nTemps); % these are equidistant
temperatures = phi(lambda);
options.llh.prior.flag = false; % compute only likelihoood
f = zeros(numel(temperatures),1);  % transformed expected log deviance
timestart1 = tic;
for l = [2:17]
    sum_ll = 0;
    thin = ceil(numel(parameters_l{l}.S.burnin+1:...
       find(~isnan(parameters_l{l}.S.logPost(:,1)),1,'last'))/1000);
    usedparams = 0;
    for j=parameters_l{l}.S.burnin+1:thin:find(~isnan(parameters_l{l}.S.logPost(:,1)),1,'last')
        ll = logLikelihood(parameters_l{l}.S.par(:,j),M,D,options.llh,conditions);
        if ~isfinite(ll)
            disp('error')
        else
            sum_ll= sum_ll+ll;
            usedparams =usedparams +1;
        end
    end
    f(l) = phiprime(lambda(l))*sum_ll/usedparams;
end
f(1) = 0;
h = lambda(3) - lambda(1);
Q = (h/6)*(f(1) + f(end) + sum(2*f(3:2:end-2)) + 4*sum(f(2:2:end-1)));
save(['./results/results_BFchains_' num2str(icomb) ],'f','h','Q')

logmarg = Q;

