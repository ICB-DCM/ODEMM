function [] = model_SP_all()
% This function generates model and estimates parameters for the model accounting for
% means and variances (using Sigma Points (SP). \n
% All kinetic parameters are assumed to vary between individual cells: \n
% k1: inter- and intra-subpopulation variable \n
% k2: cell-to-cell variable \n
% k3: cell-to-cell variable \n

load('./data/conversionprocess_data')
parameters.name = {'m(log(k_{1,s1}))','m(log(k_{1,s2}))','m(log(k_2))',...
    'm(log(k_3))','log(\sigma_{k_1}^2)','log(\sigma_{k2}^2)',...
    'log(\sigma_{k_3}^2)','log_{10}(\sigma_{noise}))',......
    'w_1'};

parameters.number = length(parameters.name);
parameters.min = [-3*ones(8,1);0];
parameters.max = [ 3*ones(4,1);2;2;2;2;1];
%%
M.name = 'CR_SP_all';
M.n_subpop = 2; %number of subpopulations
M.model = @(T,theta,u) simulate_SP_CR_all(T,theta,[]);
e = 1; % experiment index
for s = 1:M.n_subpop
    M.mean_ind{s,e} = [1]; %index of mean in output
    M.var_ind{s,e} = [2]; %index of variances in output (empty if RREs used)
end
M.sim_type = 'HO';

% definition of symbolic parameters required for generateODEMM
xi = sym(zeros(parameters.number,1));
for i = 1:parameters.number
    xi(i) = sym(['xi_' num2str(i,'%d')]);
end
% definition of input/differences
u = sym(zeros(1,1)); % 1 for definition of stimulation, 1 for subpop, 1 for experiment
for i = 1:2
    u(i) = sym(['u_' num2str(i,'%d')]);
end
r=1;
n_siminput = size(D(e).u,1); % number of input for simulation of inidividual subpopulation

% definition of parameters for simulation
M.sym.theta = [log(10)*(u(n_siminput+1)*xi(1)+(1-u(n_siminput+1))*xi(2));... %k1
    log(10)*xi(3);...
    log(10)*xi(4);...
    log(10)*2*xi(5);...
    log(10)*2*xi(6);...
    log(10)*2*xi(7)];

M.u{1,e} = [1]; % subpop
M.u{2,e} = [0];

% distribution assumption
for s = 1:M.n_subpop
    M.distribution{s,e} = 'logn_mean'; % 'norm', 'logn_median'
end
% define weights
M.sym.w{1,e} = ones(size(D(e).t(:)))*xi(9);
M.sym.w{2,e} = ones(size(D(e).t(:)))*(1-xi(9));

% define scaling, offset and measurement noise
M.sym.scaling{r,e} = sym('1'); % no scaling
M.sym.offset{r,e} = sym('0'); % no offset
M.sym.sigma_noise{r,e} = [10.^(xi(8))];

[conditions,D] = collectConditions(D,M);
options.measurement_noise = true;
options.noise_model = 'multiplicative';

generateODEMM(D,M,parameters,conditions,options);
eval(['ODEMM_' M.name]);

options.simulate_musigma = 1;
xi = (parameters.max+parameters.min)/2;
[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logLikelihood(xi,M,D,options,conditions),1e-4);
[g,g_fd_f,g_fd_b,g_fd_c]


options.MS = PestoOptions();
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('GradObj','on','display','iter',...
    'TolFun',1e-10, 'TolX',1e-10, 'MaxIter', 1000,'algorithm','interior-point');
options.MS.n_starts = 50;
options.MS.comp_type = 'sequential'; options.MS.mode = 'visual';
parameters.guess =   getParameterGuesses(parameters,@(xi) logLikelihood(xi,M,D,options,conditions),...
    options.MS.n_starts, parameters.min,parameters.max);
parameters = getMultiStarts(parameters,@(xi) logLikelihood([xi],M,D,options,conditions),options.MS);
save('./results/results_SP_all','M','D','parameters','conditions','options')
end