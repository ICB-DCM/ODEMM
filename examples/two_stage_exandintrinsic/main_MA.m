% This script generates the model incorporating the moment approximation, 
% estimates the parameters and calculates the profile likelihoods.

clear all;
close all;
clc;

%% Measurement data
load('./data/data_geneExp');

%% Definition parameters
parameters.name = {'log_{10}(k_1)','log_{10}(k_{2,1})','log_{10}(k_{2,2})',...
    'log_{10}(k_{3})',...
    'log_{10}(k_4)','log_{10}(k_5)','w'};

parameters.number = length(parameters.name);
parameters.min = [-3*ones(6,1);0];
parameters.max = [ 3*ones(6,1);1];

%% Model definition
M.name = 'geneExp_MA';
M.n_subpop = 2; %number of subpopulations
M.model = @(T,theta,u) simulate_geneExp_MA(T,theta,[]);

e = 1; % experiment index
for s= 1:M.n_subpop
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

n_siminput = size(D(e).u,1); % number of input for simulation of inidividual subpopulation
M.sym.theta = [10.^xi(1);... %k1
    u(n_siminput+1)*10.^xi(2)+(1-u(n_siminput+1))*10.^xi(3);... %k2
    10.^xi(4);... %k3
    10.^xi(5);...%k4
    10.^xi(6)];%k5

M.u{1,e} = [1]; % subpop, condition, experiment
M.u{2,e} = [0];

% distribution assumption
for s = 1:M.n_subpop
    M.distribution{s,e} = 'norm'; 
end
% initialize weights
M.sym.w{1,e} = ones(size(D(e).t(:)))*xi(7);
M.sym.w{2,e} = ones(size(D(e).t(:)))*(1-xi(7));

r=1;
M.sym.scaling{r,e} = [sym('1')];
M.sym.offset{r,e} = [sym('0')];

[conditions,D] = collectConditions(D,M);

generateODEMM(D,M,parameters,conditions)
eval(['ODEMM_' M.name]);

%% Multi-start optimization
options.MS = PestoOptions();
options.MS.localOptimizerOptions = optimset('GradObj','on',...
    'display','iter','TolFun',1e-10, ...
    'TolX',1e-10, 'MaxIter', 1000,'algorithm','interior-point');
options.MS.n_starts = 100;
options.MS.comp_type = 'sequential'; options.MS.mode = 'visual';
parameters = getMultiStarts(parameters,@(xi) logLikelihood(xi,M,D,options,conditions),options.MS);
parameters.MS.BIC = -2*parameters.MS.logPost+ log(sum(sum(~isnan(D(1).y))))*parameters.number;
save ./results/geneExp_MA parameters D M options conditions

%% Plot Fit
load_plot_settings
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.data.col{1} = color.data;
options_plot.boundaries(1).y_min = 0;
options_plot.boundaries(1).y_max = 2;
xi = parameters.MS.par(:,1);
options_plot.model.col{1} = color.SP;
fh = plotODEMM(D,M,xi,options_plot);

%% Profile likelihood calculation
options.MS.localOptimizerOptions = optimset('GradObj','on',...
    'display','off','MaxIter',50,...
    'algorithm','trust-region-reflective');
parameters = getParameterProfiles(parameters,@(xi) ...
    logLikelihood(xi,M,D,options,conditions),options_PL);
save ./results/geneExp_MA parameters D M options conditions

plotParameterProfiles(parameters)

