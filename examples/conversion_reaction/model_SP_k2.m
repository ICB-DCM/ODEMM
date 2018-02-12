function [] = model_SP_k2()
% This function generates the model and estimates the parameters in case of 
% accounting for  means and variances (using Sigma Points (SP)).
% The parameter are assumed to be: \n
% k1: subpopulation variable  \n
% k2: cell-to-cell variable  \n
% k3: homogeneous

load('./project/data/conversionprocess_data')

parameters.name = {'log_{10}(k_{1,1})','log_{10}(k_{1,2})','log_{10}(k_2)',...
                   'log_{10}m(log(k_3))','log_{10}(\sigma_{k_3}))','log_{10}(\sigma_{noise}))',......
                    'w'};               
parameters.number = length(parameters.name);
parameters.min = [-3*ones(4,1);-3;-3;0];
parameters.max = [ 3*ones(4,1);2;2;1];

%% Model definition and generation
M.name = 'CR_SP_k2'; 
M.n_subpop = 2; %number of subpopulations
M.model = @(T,theta,u) simulate_SP_CR_k2(T,theta,[]);
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

u = sym(zeros(1,1)); % 1 for definition of stimulation, 1 for subpop, 1 for experiment
for i = 1:2
   u(i) = sym(['u_' num2str(i,'%d')]);
end
r=1;
n_siminput = size(D(e).u,1); % number of input for simulation of inidividual subpopulation
M.sym.theta = [log(10)*(u(n_siminput+1)*xi(1)+(1-u(n_siminput+1))*xi(2));... %k1             
               log(10)*xi(3);...
               log(10)*xi(4);...
               log(10)*2*xi(5)];
                     
M.u{1,e} = [1]; % subpop
M.u{2,e} = [0];
      
% distribution assumption
for s = 1:M.n_subpop
    M.distribution{s,e} = 'logn_mean'; % 'norm', 'logn_median'
end
% define weights
M.sym.w{1,e} = ones(size(D(e).t(:)))*xi(7);
M.sym.w{2,e} = ones(size(D(e).t(:)))*(1-xi(7));
% define scaling, offset and measurement noise
M.sym.scaling{r,e} = sym('1');
M.sym.offset{r,e} = sym('0');
M.sym.sigma_noise{r,e} = [10.^(xi(6))];
options.measurement_noise = true;
options.noise_model = 'multiplicative';
%% generate model

[conditions,D] = collectConditions(D,M);
generateODEMM(D,M,parameters,conditions,options);
eval(['ODEMM_' M.name]);

options.simulate_musigma = 1;
xi = (parameters.max+parameters.min)/2;
[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logLikelihood(xi,M,D,options,conditions),1e-4);
[g,g_fd_f,g_fd_b,g_fd_c]

%% Multi-start optimization
options.MS = PestoOptions();
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('GradObj','on','display','iter',...
    'TolFun',1e-10, 'TolX',1e-10, 'MaxIter', 1000,'algorithm','interior-point');
options.MS.n_starts = 50;
options.MS.comp_type = 'sequential'; options.MS.mode = 'visual';

parameters.guess =   getParameterGuesses(parameters,@(xi) logLikelihood(xi,M,D,options,conditions),...
    options.MS.n_starts, parameters.min,parameters.max);
parameters = getMultiStarts(parameters,@(xi) logLikelihood([xi],M,D,options,conditions),options.MS);
parameters.MS.BIC = -2*parameters.MS.logPost+ log(numel(D(1).t)*1000)*parameters.number;
save('./project/results/results_SP_k2','M','D','parameters','conditions','options')