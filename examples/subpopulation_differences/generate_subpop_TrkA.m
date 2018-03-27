% Generation of final model including differences in TrkA levels.

clear all
close all
clear all

load('./data/data_PDL');

parameters.max = [6*ones(8,1);... %kinetic
    6*ones(11,1);...%scaling
    6*ones(9,1); % noise
    1]; % weight
parameters.min = [-6*ones(8,1);... %kinetic
    -6*ones(11,1);...%scaling
    -6*ones(9,1); % noise
    0]; %

parameters.name = {'log_{10}(k_1)',...
'log_{10}(k_2)',...
'log_{10}(k_4)',...
'log_{10}(k_5)',...
'log_{10}(\beta_{k_3[TrkA]_{0,1}})',...
'log_{10}(\beta_{k_3[TrkA]_{0,2}})',...
'log_{10}(\beta_{c_P^{1}[Erk]_{0}})',...
'log_{10}(\sigma_{T,1})',...
'log_{10}(\sigma_{T,2})',...
'log_{10}(\sigma_{E})',...
'log_{10}(\sigma_{TE})',...
'log_{10}(o_P^{1})',...
'log_{10}(c_P^{2})',...
'log_{10}(o_P^{2})',...
'log_{10}(c_P^{3})',...
'log_{10}(o_P^{3})',...
'log_{10}(c_P^{4})',...
'log_{10}(o_P^{4})',...
'log_{10}(c_T/k_3)',...
'log_{10}(o_T)',...
'log_{10}(c_E/c_P^{1})',...
'log_{10}(o_E)',...
'log_{10}(\sigma_{P,noise}^{1})',...
'log_{10}(\sigma_{P,noise}^{2})',...
'log_{10}(\sigma_{P,noise}^{3})',...
'log_{10}(\sigma_{P,noise}^{4})',...
'log_{10}(\sigma_{T,noise})',...
'log_{10}(\sigma_{E,noise})',...
'w_1'};

M.n_subpop = 2; %number of subpopulations
M.model = @(T,theta,u)simulate_SigmaPoints_sPsET_loglog_corr(T,theta,u(1));
%
% P E T PP PE PT EE ET TT
for s = 1:M.n_subpop
    for e=1:2
        M.mean_ind{s,e} = [1]; %index of mean in output
        M.var_ind{s,e} = [4]; %index of variances in output (empty if RREs used)
    end
    e = 3;
    M.mean_ind{s,e} = [3,1]; %TrkA/pErk
    M.var_ind{s,e} = [9,6,4];
    e = 4;
    M.mean_ind{s,e} = [2,1]; %Erk/pErk
    M.var_ind{s,e} = [7,5,4];
    
end

M.sim_type = 'HO';

options.replicates = false;
options.write_parameters = true;
options.measurement_noise = true;
options.noise_model = 'multiplicative';
options.penalize_scaling = false;
comb = [0,0,0,0,1,1];

parameters.number = length(parameters.name);

% definition of symbolic parameters required for generateODEMM
xi = sym(zeros(parameters.number,1));
for i = 1:parameters.number
    xi(i) = sym(['xi_' num2str(i,'%d')]);
end
u = sym(zeros(2,1)); % 1 for definition of stimulation, 1 for subpop,
for i = 1:7
    u(i) = sym(['u_' num2str(i,'%d')]);
end
g = size(D(e).u,1); % number of input for simulation of inidividual subpopulation

for e=[1:4]
    M.u{1,e} = [0]; % subpop
    M.u{2,e} = [1];
end
M.sym.theta = [log(10)*(xi(1));... %k1
    log(10)*(xi(2));... %k2
    log(10)*(xi(3));... %k4
    log(10)*(xi(4));... %k5
    log(10)*(u(g+1)*xi(5) + (1-u(g+1))*xi(6));... %k3TrkA
    log(10)*(xi(7));...  %sPmErk
    log(log(10.^(2*(u(g+1)*xi(8) + (1-u(g+1))*xi(9)))+1));... % log(sigma^2)
    log(10)*xi(11);...
    log(log(10.^(2*(xi(10)))+1))];
 
r=1;
e=1;
M.sym.scaling{r,e} = sym('1');
M.sym.offset{r,e} = 10.^(xi(12));
M.sym.sigma_noise{e} = 10.^(xi(23));

e=2;
M.sym.scaling{r,e} = 10.^(xi(13));
M.sym.offset{r,e} = 10.^(xi(14));
M.sym.sigma_noise{e} = 10.^(xi(24));

e=3;
M.sym.scaling{r,e} = [10.^xi(19);10.^xi(15)];
M.sym.offset{r,e} = [10.^xi(20);10.^xi(16)];
M.sym.sigma_noise{e} = [10.^xi(27);10.^xi(25)];

e=4;
M.sym.scaling{r,e} = [10.^xi(21);10.^xi(17)];
M.sym.offset{r,e} = [10.^xi(22);10.^xi(18)];
M.sym.sigma_noise{e} = [10.^xi(28);10.^xi(26)];

for e = 1:4
    for s = 1:2
        M.distribution{s,e} = 'logn_mean'; 
    end
    M.sym.w{1,e} = ones(size(D(e).t(:)))*xi(29);
    M.sym.w{2,e} = ones(size(D(e).t(:)))*(1-xi(29));
end
M.name = ['NGF_subpop_TrkA'];
options.simulate_musigma = true;
[conditions,D] = collectConditions(D,M);
generateODEMM(D,M,parameters,conditions,options);
eval(['ODEMM_' M.name])
