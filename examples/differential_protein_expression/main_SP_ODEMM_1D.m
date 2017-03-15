% This script generates the model for the 1-dimensional measurements of the
% differential gene expression example using the sigma point approximation.
% The parameters are estimated and also the profile likelihoods are
% calculated using the toolbox PESTO.

clear all
close all
clc

%% load data
load('./project/data/oneStage_data_1D')

%% Parameter definition
parameters.name = {'log_{10}(\lambda_0)','log_{10}(\lambda_{A,1})','log_{10}(\lambda_{A,2})',...
    'log_{10}(\lambda_{B,1})','log_{10}(\lambda_{B,2})','log_{10}(m\gamma)',...
    'log_{10}(\Sigma_{\gamma})','log_{10}(\sigma_{noise})'...
    'w_1'};
parameters.number = length(parameters.name);
parameters.min = [-3*ones(6,1);-3;-3;0];
parameters.max = [ 3*ones(6,1);1;1;1];

%% Model definition and generation
M.name = 'oneStage_SP_1D';
M.n_subpop = 2; % number of subpopulations
M.model = @(T,theta,u) simulate_onestage_SP(T,theta,[]);
for s = 1:M.n_subpop
    e=1; % first data set
    M.mean_ind{s,e} = 1; %index of mean in output
    M.var_ind{s,e} = [3]; %index of variances in output (empty if RREs used)
    e=2; % second dataset
    M.mean_ind{s,e} = 2; %index of mean in output
    M.var_ind{s,e} = [5]; %index of variances in
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
M.sym.theta = [xi(1)*log(10);...
    log(10)*(u(n_siminput+1)*xi(2)+(1-u(n_siminput+1))*xi(3));...
    log(10)*(u(n_siminput+1)*xi(4)+(1-u(n_siminput+1))*xi(5));... %k4
    log(10)*(xi(6));...
    log(10)*2*xi(7)];

for e = 1:numel(D)
    for s = 1:M.n_subpop
        % distribution assumption
        M.distribution{s,e} = 'logn_mean'; % 'norm', 'logn_median'
    end
    % definition weights/difference indicatos
    s=1;
    M.sym.w{s,e} = ones(size(D(e).t(:)))*xi(9);
    M.u{s,e} = [1]; 

    s=2;
    M.sym.w{s,e} = ones(size(D(e).t(:)))*(1-xi(9));
    M.u{s,e} = [0];

    % definition scaling, offset and measurement noise
    M.sym.scaling{r,e} = [sym('1')];
    M.sym.offset{r,e} = [sym('0')];
    M.sym.sigma_noise{r,e} = [10.^(xi(8))];
end

options.measurement_noise = true;
options.noise_model = 'multiplicative';
options.simulate_musigma = 1;

[conditions,D] = collectConditions(D,M);

generateODEMM(D,M,parameters,conditions,options);
eval(['ODEMM_' M.name]);

xi = (parameters.max+parameters.min)/2;
[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logLikelihood(xi,M,D,options,conditions),1e-4);
[g,g_fd_f,g_fd_b,g_fd_c]

%% Multi-start optimization
options.MS.fmincon = optimset('GradObj','on','display','iter','TolFun',1e-10, 'TolX',1e-10, 'MaxIter', 1000,'algorithm','interior-point');
options.MS.n_starts = 100;
options.MS.comp_type = 'sequential'; options.MS.mode = 'visual';
parameters = getMultiStarts(parameters,@(xi) logLikelihood([xi],M,D,options,conditions),options.MS);
save('./project/results/results_diffgeneexp_1D','M','D','options',...
     'parameters','parameters','conditions')
%% Profile likelihood calculation for both minima (global and local)
options.PL.fmincon = optimset('GradObj','on','display','off','MaxIter',100,'algorithm','interior-point');
options.PL.parameter_index = 1:parameters.number;
options.PL.P.min = max(-6,parameters.min);
options.PL.P.max = min( 6,parameters.max);
options.PL.P_next_step.min = 1e-5;
options.PL.MAP_index = 1;
parameters = getParameterProfiles(parameters,@(xi) logLikelihood(xi,M,D,options,conditions),options.PL);
options.PL.parameter_index = [2,3,4,5];
options.PL.MAP_index = 43;
parameters_2ndmode = getParameterProfiles(parameters,@(xi) logLikelihood(xi,M,D,options,conditions),options.PL);
save('./project/results/results_diffgeneexp_1D','M','D','options',...
    'parameters_2ndmode','parameters','conditions')

%% Profile likelihood calculation using Profile Integration
% options.PL.solver = struct('gamma', 0, ...
%     'type', 'ode45', ...'ode45', 'ode15s', 'ode15sDAE', 'ode113','CVODE'
%     'algorithm', 'Adams', ...
%     'nonlinSolver', 'Newton', ...
%     'linSolver', 'Dense', ...
%     'eps', 1e-7, ...
%     'minCond', 1e-7, ...
%     'gradient', true, ...
%     'hessian', 'fd', ...
%     'MaxStep', 1, ...
%     'MinStep', 1e-3, ...
%     'MaxNumSteps', 1e4, ...
%     'RelTol', 1e-3, ...
%     'AbsTol', 1e-4, ...
%     'hessianStep', 1e-5, ...
%     'gradientStep', 1e-5 ...
%     );
% options.PL.mode     = 'visual';
% options.PL.fh = [];
% options.PL.save = 0;
% options.PL.plot_options = [];
% parameters = getProfileIntegration(parameters, @(xi) logLikelihood(xi,M,D,options,conditions), options.PL);