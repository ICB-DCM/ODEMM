clear all
close all
clc

%% Data generation 
k1 = 10; % basal mA
k2 = 10; % mA stimulus subpop 1
k2_2 = 20; % mA stimulus subpop 2
k3 = 1; % degradation mA
k4 = 5; % production A
k5 = 0.1; % degradation A
w = 0.3;

load data/data_geneExp1D

parameters.name = {'log_{10}(k_{1})',...
                   'log_{10}(k_{2,1})',...
                   'log_{10}(k_{2,2})',...
                   'log_{10}(k_3)',...
                   'log_{10}(k_4)',...
                   'log_{10}(k_5)',...
                   'w'};
parameters.number = length(parameters.name);
M.n_subpop = 2; %number of subpopulations
M.model = @(T,theta,u) simulate_geneExp1D(T,theta,[]); 
e=1; % experiment index (only one experiment)
for s = 1:2
    M.mean_ind{s,e} = [1]; %index of mean in output
    M.var_ind{s,e} = [2]; %index of variances in output (empty if RREs used)
end
M.sim_type = 'HO';
xi = sym(zeros(parameters.number,1));
for i = 1:parameters.number
    xi(i) = sym(['xi_' num2str(i,'%d')]);
end
n_siminput = size(D(1).u,1); % number of input for simulation of inidividual subpopulation
u = sym(zeros(n_siminput+1,1));
for i = 1:n_siminput+1
    u(i) = sym(['u_' num2str(i,'%d')]);
end
M.sym.theta = [10.^xi(1);u(n_siminput+1)*10.^xi(2)+(1-u(n_siminput+1))*10.^xi(3);10.^xi(4);10.^xi(5);10.^xi(6)];
M.u{1,e} = [1]; % subpop 1
M.u{2,e} = [0]; % subpop 2
% initialize weights
M.sym.w{1,e} = xi(7)*ones(size(D(e).t(:)));
M.sym.w{2,e} = (1-xi(7))*ones(size(D(e).t(:)));
r=1; % only one replicate
M.sym.scaling{r,e} = sym('1');
M.sym.offset{r,e} = sym('0');

% distribution assumption
for s = 1:M.n_subpop
    M.distribution{s,e} = 'norm';
end

[conditions,D] = collectConditions(D,M);
options.dimension = 'univariate'; % 1D-measurements

parameters.min = [-3*ones(6,1);0];
parameters.max = [ 3*ones(6,1);1];

M.name = [M.distribution{1,1} '_geneExp1D'];

generateODEMM_extended(D,M,parameters,conditions,options);
eval(['ODEMM_' M.name]);

if strcmp(M.distribution{s,e},'students_t')
    parameters.max(end) = 8;
end
%%
xi_true = [k1,k2,k2_2,k3,k4,k5,w];
xi_true(1:6) = log10(xi_true(1:6));
options.use_robust = 1;
xi = rand(parameters.number,1);
ll = logLikelihood_extend(xi,M,D,options);
[ll,dll] =  logLikelihood_extend(xi,M,D,options);
[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logLikelihood_extend(xi,M,D,options),1e-6);
[g,g_fd_f,g_fd_b,g_fd_c]

%%
options.MS = PestoOptions();
options.MS.n_starts = 10;
warning off
parameters.guess =   getParameterGuesses(parameters,@(xi) logLikelihood_extend(xi,M,D,options,conditions),...
    options.MS.n_starts, parameters.min,parameters.max);
parameters = getMultiStarts(parameters,@(xi) logLikelihood_extend([xi],M,D,options,conditions),options.MS);
xi = parameters.MS.par(:,1)
plotODEMM(D,M,xi)

%%
xi_true(end+1) = log10(3);
plotODEMM(D,M,xi_true)
