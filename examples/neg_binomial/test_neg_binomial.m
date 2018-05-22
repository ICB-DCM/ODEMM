clear all
close all
clc

%% Data generation
t = [0,0.1,0.5,1,2]; %hours
xi_true = [2,2,0,0.5];

theta = xi_true(1:3); % log
w = 1; % lin
n_data = 1e3;
D(1).y = nan(1,numel(t),n_data,1);
D(1).t = t;
D(1).n_dim = 1;
D(1).measurand = 'conc. of B';
D(1).u = 1;
D(1).name = 'conversion reaction';
rho = xi_true(4);
for i = 1:w*n_data
    sol = simulate_onestage_1D(t,theta,[]); 
    for k = 1:numel(t)
        D(1).y(1,k,i,:) = nbinrnd((1-rho)/rho*sol.y(k),rho);
    end
end
% for i = w*n_data+1:n_data
%     i
%     theta_s2(3) = (m_k3+randn(1)*sigma_k3);
%     sol = simulate_CR_log(t,theta_s2,[]);
%     D(1).y(1,:,i,:) = exp(sol.y'+randn(size(sol.y'))*sigma_noise);
%     trajectories(:,i) = exp(sol.y);
% end
options_plot.data.bins = 20;
plotODEMM(D)



%% save data
save data/negbin_data D xi_true

%%

parameters.name = {'log_{10}(k_{1})',...
    'log_{10}(k_{2})','log_{10}(k_3)'};
parameters.number = length(parameters.name);

M.n_subpop = 1; %number of subpopulations
M.model = @(T,theta,u) simulate_onestage_1D(T,theta,[]); 

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
M.sym.theta = [xi(1);xi(2);xi(3)];
M.u{1,e} = [1]; % subpop

% distribution assumption
for s = 1:M.n_subpop
    M.distribution{s,e} = 'neg_binomial';
end
% initialize weights
M.sym.w{1,e} = sym(ones(size(D(e).t(:))));
r=1; % only one replicate
M.sym.scaling{r,e} = sym('1');
M.sym.offset{r,e} = sym('0');
[conditions,D] = collectConditions(D,M);
options.dimension = 'univariate'; % 1D-measurements (for getRREsigmas)

parameters.min = [-3*ones(3,1)];
parameters.max = [ 3*ones(3,1)];

options.rhos = 'time-dependent';
[parameters,conditions,D] = getRRErhos(parameters,conditions,options,D,M);
M.name = 'test_neg_binomial';
generateODEMM_extended(D,M,parameters,conditions,options);
eval(['ODEMM_' M.name]);

%%
options.use_robust = 1;
xi = rand(parameters.number,1);
ll = logLikelihood_extend(xi,M,D,options);
[ll,dll] =  logLikelihood_extend(xi,M,D,options);
[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logLikelihood_extend(xi,M,D,options),1e-6);
[g,g_fd_f,g_fd_b,g_fd_c]

plotODEMM(D,M,xi_true)
