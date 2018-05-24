% This script generates the data for the two stage gene expression example
% using the Stochastic Simulation Algorithm (SSA) implemented in the
% toolbox CERENA.

clear all
clc
close all

modelDefFileName = 'modelDef_geneExp';
eval(modelDefFileName);
System = completeSystem(System);
System = completeSystemSSA(System);
options.mode = 'constant';

%% Data generation
k1 = 10; % basal mA/mB
k2 = 10; % mA stimulus
k2_2 = 20; % mB stimulus
k3 = 1; % deg mA/mB
k4 = 5; % A/B
k5 = 0.1; % deg A/B
w = 0.3;
t = [0,5,10,30,60];


D(1).t = t;
D(1).u = 1;
D(1).y = NaN(1,numel(t),1000);
D(1).n_dim = 1;
D(1).measurand = {'A'};
D(1).name = 'Gene expression';
n_cells = 1000;

% SSA runs for the initial conditions (without stimulation)
theta_SSA_forinit = [k1;0;k3;k4;k5];
[Xinits,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = simulate_SSA(System,...
    [1,10,300],theta_SSA_forinit,[],5000);

% Generate first population
theta_SSA = [k1;k2;k3;k4;k5];
for k = 1:numel(t)
    for j = 1:round(w*n_cells)
        System = modelDef_geneExp_changeinit(Xinits(end,:,randi([1,5000],1,1))');
        System = completeSystem(System);
        System = completeSystemSSA(System);
        [X_ssa,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = ...
            simulate_SSA(System,t(k),theta_SSA,[],1);
        if noiseflag
            
        end
        D(1).y(1,k,j) = Y_ssa;
    end
end

% Generation second population
theta_SSA = [k1;k2_2;k3;k4;k5];
for k = 1:numel(t)
    for j = round(w*n_cells)+1:n_cells
        System = modelDef_geneExp_changeinit(Xinits(end,:,randi([1,5000],1,1))');
        System = completeSystem(System);
        System = completeSystemSSA(System);
        [X_ssa,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = ...
            simulate_SSA(System,t,theta_SSA,[],1);

        D(1).y(1,k,j) = permute(Y_ssa,[1,3,2]);
    end
end


save data_geneExp_2subpop D


