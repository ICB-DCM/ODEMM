function [] = models_RRE()
% This function generates Reaction Rate Euqation model and estimates 
% the parameters for the models with varying number of parameters for the variances

load('./data/conversionprocess_data')

parameters.name = {'log_{10}(k_{1,1})',...
    'log_{10}(k_{1,2})','log_{10}(k_{2})','log_{10}(k_3)',...
    'w'};
parameters.number = length(parameters.name);

M.n_subpop = 2; %number of subpopulations
M.model = @(T,theta,u) simulate_CR(T,theta,[]);

e=1; % experiment index (only one experiment)
for s = 1:2
    M.mean_ind{s,e} = [1]; %index of mean in output
    M.var_ind{s,e} = []; %index of variances in output (empty if RREs used)
end
M.sim_type = 'RRE';

% definition of symbolic parameters required for generateODEMM
xi = sym(zeros(parameters.number,1));
for i = 1:parameters.number
    xi(i) = sym(['xi_' num2str(i,'%d')]);
end
% definition of stimulations/differences
u = sym(zeros(1,1)); % 1 for definition of stimulation, 1 for subpop, 1 for experiment
for i = 1:2
    u(i) = sym(['u_' num2str(i,'%d')]);
end
n_siminput = size(D(e).u,1); % number of input for simulation of inidividual subpopulation
M.sym.theta = [log(10)*(u(n_siminput+1)*xi(1)+(1-u(n_siminput+1))*xi(2));... %k1
    log(10)*xi(3);...
    log(10)*xi(4)];

M.u{1,e} = [1]; % subpop
M.u{2,e} = [0];

% distribution assumption
for s = 1:M.n_subpop
    M.distribution{s,e} = 'logn_mean';
end
% initialize weights
M.sym.w{1,e} = ones(size(D(e).t(:)))*xi(5);
M.sym.w{2,e} = ones(size(D(e).t(:)))*(1-xi(5));

r=1; % only one replicate
M.sym.scaling{r,e} = sym('1');
M.sym.offset{r,e} = sym('0');

[conditions,D] = collectConditions(D,M);

options.dimension = 'univariate'; % 1D-measurements (for getRREsigmas)

% options for multi-start optimization
options.MS = PestoOptions();
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('GradObj','on','display','iter','TolFun',1e-10, ...
    'TolX',1e-10, 'MaxIter', 1000,'algorithm','interior-point');
options.MS.n_starts = 50;
options.MS.comp_type = 'sequential'; 
options.MS.mode = 'visual';


model_names = {'onlyone','subpop','timedep'};
sigma_string = {'only-one','subpopulation-specific','time-dependent'};

for model_ind = 1:3
    
    parameters.name = {'log_{10}(k_{1,1})','log_{10}(k_{1,2})',...
        'log_{10}(k_{2})','log_{10}(k_3)','w'};
    parameters.number = length(parameters.name);
    parameters.min = [-3*ones(4,1);0];
    parameters.max = [ 3*ones(4,1);1];
    options.sigmas = sigma_string{model_ind};
    [parameters,conditions,D] = getRREsigmas(parameters,conditions,options,D,M);
    M.name = ['CR_RRE_' model_names{model_ind}];
    generateODEMM(D,M,parameters,conditions,options);
    eval(['ODEMM_' M.name]);
    
    %xi = (parameters.max+parameters.min)/2;
    %[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logLikelihood(xi,M,D,options,conditions),1e-4);
    %[g,g_fd_f,g_fd_b,g_fd_c]
    
    parameters.guess =   getParameterGuesses(parameters,@(xi) logLikelihood(xi,M,D,options,conditions),...
        options.MS.n_starts, parameters.min,parameters.max);
    
    parameters = getMultiStarts(parameters,@(xi) logLikelihood([xi],M,D,options,conditions),options.MS);
    parameters.MS.BIC = -2*parameters.MS.logPost(1)+log(1000*numel(D(1).t))*parameters.number;
      save(['./results/results_RRE_' model_names{model_ind}],'M','D','parameters','conditions','options')
    clear parameters 
end