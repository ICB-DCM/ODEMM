% This script generates the data for the gene expression example
% using the Stochastic Simulation Algorithm (SSA) implemented in the
% toolbox CERENA.
clear all
clc
close all

warning off
modelDefFileName = 'modelDef_geneExp1D';
modelName = 'modelDef_geneExp1D';

method = 'MEC_2_LD_2_a';
%System = genmexp(modelName,'modelDef_geneExp1D_foramici',method);
%amiwrap('geneExp1D',[method,'_',modelName,'_syms'])
 
eval(['System = ' modelDefFileName '([])']);
System = completeSystem(System);
System = completeSystemSSA(System);
 
options.mode = 'constant';

%% Data generation 
k1 = 10; % basal mA
k2 = 10; % mA stimulus subpop 1
k2_2 = 20; % mA stimulus subpop 2
k3 = 1; % degradation mA
k4 = 5; % production A
k5 = 0.1; % degradation A

w = 0.3;

n_cells=100;
%%
theta_SSA_forinit = [k1;0;k3;k4;k5];
% SSA runs for the initial conditions (without stimulation)
[Xinits,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = simulate_SSA(System,[1,10,300],...
    theta_SSA_forinit,[],n_cells);

t = [0,5,10,30,60];
% Generate first population
D(1).t = t;
D(1).u = 1;
D(1).y = NaN(1,numel(t),n_cells);
theta_SSA = [k1;k2;k3;k4;k5];
for k = 1:round(n_cells*w)
    System = modelDef_geneExp1D(Xinits(end,:,k)');
    System = completeSystem(System);
    System = completeSystemSSA(System);
    [X_ssa,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = simulate_SSA(System,t,theta_SSA,[],1);
    D(1).y(1,:,k) = permute(Y_ssa,[1,3,2]); 
end

% Generation second population
theta_SSA = [k1;k2_2;k3;k4;k5];
for k = (round(n_cells*w)+1):n_cells
    System = modelDef_geneExp1D(Xinits(end,:,k)');
    System = completeSystem(System);
    System = completeSystemSSA(System);
    [X_ssa,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = simulate_SSA(System,t,theta_SSA,[],1);
    D(1).y(1,:,k) = permute(Y_ssa,[1,3,2]); 
end
D(1).n_dim = 1;
D(1).measurand = {'A'};
D(1).name = 'Gene expression';

plotODEMM(D)
%%
save data_geneExp1D D

sol = simulate_geneExp1D(D.t,theta_SSA);
%%
figure
theta_SSA = [k1;k2;k3;k4;k5];
sol = simulate_geneExp1D(D.t,theta_SSA);
plot(sol.t,sol.y(:,1),'r-'); hold on;
plot(D.t, squeeze(mean(D(1).y(1,:,1:30),3)),'b--');
theta_SSA = [k1;k2_2;k3;k4;k5];
sol = simulate_geneExp1D(D.t,theta_SSA);
plot(sol.t,sol.y(:,1),'r-'); hold on;
plot(D.t, squeeze(mean(D(1).y(1,:,31:end),3)),'b--');