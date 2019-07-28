% This script generates the (outlier-free) data for the simulation study.

clear all
close all
clc

% Models included in the manuscript number 1, 3 and 4
modelnames = {'conversionReaction',... % 1D data
    'diffProteinExpression',... % 2D data
    'twoStageGeneExpression',... % 1D data
    'diffProteinExpression'}; % 1D data

load_simStudy_settings

% Conversion 
M{1}.theta = @(xi,u) [u(1)*10.^(u(2)*xi(1)+(1-u(2))*xi(2));...
    10.^xi(3);...
    10.^xi(4)];
M{1}.w = @(xi,u) u(2)*xi(5)+(1-u(2))*xi(5);

% Birth-death (diffProteinExpression)
M{2}.theta = @(xi,u) [10.^xi(1);...
    u(1)*10.^(u(2)*xi(2)+(1-u(2))*xi(3));...
    u(1)*10.^(u(2)*xi(4)+(1-u(2))*xi(5));
    10.^xi(6)];
M{2}.w = @(xi,u) u(2)*xi(7)+(1-u(2))*xi(7);

% Two-stage
M{3}.theta = @(xi,u) [10.^xi(1);
    u(1)*10.^(u(2)*xi(2)+(1-u(2))*xi(3));...
    10.^xi(4);...
    10.^xi(5);
    10.^xi(6)];
M{3}.w = @(xi,u) u(2)*xi(7)+(1-u(2))*xi(7);

% Birth-death (diffProteinExpression)
M{4}.theta = @(xi,u) [10.^xi(1);...
    u(1)*10.^(u(2)*xi(2)+(1-u(2))*xi(3));
    10.^xi(4)];
M{4}.w = @(xi,u) u(2)*xi(5)+(1-u(2))*xi(5);

%% Simulate means and variances
m = 4; % example for model number 4
t = tps{1};
tsim = linspace(t(1),t(end));
xi = parametersSets{m}(3,:);
if m == 1
    sol1 = simulate_model1_conversionReaction(tsim,M{1}.theta(xi,[1,0]));
    sol2 = simulate_model1_conversionReaction(tsim,M{1}.theta(xi,[1,1]));
    plot(sol1.t, sol1.y(:,1),'r'); hold on;
    plot(sol1.t, sol1.y(:,1)+sqrt(sol1.y(:,2)),'r--'); hold on;
    plot(sol1.t, sol1.y(:,1)-sqrt(sol1.y(:,2)),'r--'); hold on;
    plot(sol2.t, sol2.y(:,1),'b'); hold on;
    plot(sol2.t, sol2.y(:,1)+sqrt(sol2.y(:,2)),'b--'); hold on;
    plot(sol2.t, sol2.y(:,1)-sqrt(sol2.y(:,2)),'b--'); hold on;
elseif m == 2
    sol1 = simulate_model2_diffProteinExpression(tsim,M{2}.theta(xi,[1,0]));
    sol2 = simulate_model2_diffProteinExpression(tsim,M{2}.theta(xi,[1,1]));
    subplot(1,2,1)
    plot(sol1.t, sol1.y(:,1),'r'); hold on;
    plot(sol1.t, sol1.y(:,1)+sqrt(sol1.y(:,3)),'r--'); hold on;
    plot(sol1.t, sol1.y(:,1)-sqrt(sol1.y(:,3)),'r--'); hold on;
    plot(sol2.t, sol2.y(:,1),'b'); hold on;
    plot(sol2.t, sol2.y(:,1)+sqrt(sol2.y(:,3)),'b--'); hold on;
    plot(sol2.t, sol2.y(:,1)-sqrt(sol2.y(:,3)),'b--'); hold on;
    subplot(1,2,2)
    plot(sol1.t, sol1.y(:,2),'r'); hold on;
    plot(sol1.t, sol1.y(:,2)+sqrt(sol1.y(:,5)),'r--'); hold on;
    plot(sol1.t, sol1.y(:,2)-sqrt(sol1.y(:,5)),'r--'); hold on;
    plot(sol2.t, sol2.y(:,2),'b'); hold on;
    plot(sol2.t, sol2.y(:,2)+sqrt(sol2.y(:,5)),'b--'); hold on;
    plot(sol2.t, sol2.y(:,2)-sqrt(sol2.y(:,5)),'b--'); hold on;
elseif m == 3
    sol1 = simulate_model3_twoStageGeneExpression(tsim,M{m}.theta(xi,[1,0]));
    sol2 = simulate_model3_twoStageGeneExpression(tsim,M{m}.theta(xi,[1,1]));
    plot(sol1.t, sol1.y(:,1),'r'); hold on;
    plot(sol1.t, sol1.y(:,1)+sqrt(sol1.y(:,2)),'r--'); hold on;
    plot(sol1.t, sol1.y(:,1)-sqrt(sol1.y(:,2)),'r--'); hold on;
    plot(sol2.t, sol2.y(:,1),'b'); hold on;
    plot(sol2.t, sol2.y(:,1)+sqrt(sol2.y(:,2)),'b--'); hold on;
    plot(sol2.t, sol2.y(:,1)-sqrt(sol2.y(:,2)),'b--'); hold on;
elseif m == 4
    sol1 = simulate_model4_diffProteinExpression(tsim,M{m}.theta(xi,[1,0]));
    sol2 = simulate_model4_diffProteinExpression(tsim,M{m}.theta(xi,[1,1]));
    plot(sol1.t, sol1.y(:,1),'r'); hold on;
    plot(sol1.t, sol1.y(:,1)+sqrt(sol1.y(:,2)),'r--'); hold on;
    plot(sol1.t, sol1.y(:,1)-sqrt(sol1.y(:,2)),'r--'); hold on;
    plot(sol2.t, sol2.y(:,1),'b'); hold on;
    plot(sol2.t, sol2.y(:,1)+sqrt(sol2.y(:,2)),'b--'); hold on;
    plot(sol2.t, sol2.y(:,1)-sqrt(sol2.y(:,2)),'b--'); hold on;
end
xlabel('t')
ylabel('count')

%% Generate data with SSA
for m = [1,3,4]
    for it = 1:size(tps,2)
        t = tps{it};
        for ic = 1:4
            for set = 1:3
                try % if it already exists, do not generate it again
                    load(['./data/data_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                        'cells_' num2str(length(t)) 'tps_' num2str(set) 'paramsetD'],'D')
                catch
                    clear D
                    xi = parametersSets{m}(set,:);
                    modelDefFileName = ['modelDef_model' num2str(m) '_' modelnames{m}];
                    modelName = ['model' num2str(m) '_' modelnames{m}];
                    eval(['System = ' modelDefFileName '([]);']);
                    System = completeSystem(System);
                    System = completeSystemSSA(System);
                    theta_SSA_forinit = M{m}.theta(xi,[0,0]);
                    
                    % SSA runs for the initial conditions (without stimulation)
                    [Xinits,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = simulate_SSA(System,[1,10,300],...
                        theta_SSA_forinit,[],n_cells(ic));
                    
                    w = M{m}.w(xi,[0,0]);
                    theta_SSA = M{m}.theta(xi,[1,0]);
                    for k = 1:round(n_cells(ic)*w)
                        k
                        eval(['System = modelDef_model' num2str(m) '_' modelnames{m} '(Xinits(end,:,k)'');']);
                        System = completeSystem(System);
                        System = completeSystemSSA(System);
                        [X_ssa,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = simulate_SSA(System,t,theta_SSA,[],1);
                        if m == 2
                            D(1).y(1,:,k,:) = permute(Y_ssa,[1,3,2]);
                        else
                            D(1).y(1,:,k) = permute(Y_ssa,[1,3,2]);
                        end
                    end
                    
                    % Generation second population
                    theta_SSA = M{m}.theta(xi,[1,1]);
                    for k = (round(n_cells(ic)*w)+1):n_cells(ic)
                        k
                        eval(['System = modelDef_model' num2str(m) '_' modelnames{m} '(Xinits(end,:,k)'');']);
                        System = completeSystem(System);
                        System = completeSystemSSA(System);
                        [X_ssa,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = simulate_SSA(System,t,theta_SSA,[],1);
                        if m == 2
                            D(1).y(1,:,k,:) = permute(Y_ssa,[1,3,2]);
                        else
                            D(1).y(1,:,k) = permute(Y_ssa,[1,3,2]);
                        end
                    end
                    D(1).distribution = 'SSA';
                    D(1).xi_true = xi;
                    D(1).t = t;
                    D(1).u = 1;
                    D(1).name = modelnames{m};
                    if m == 2
                        D(1).n_dim = 2;
                        D(1).measurand = {'A','B'};
                    else
                        D(1).n_dim = 1;
                        D(1).measurand = {'A'};
                    end
                    save(['./data/data_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                        'cells_' num2str(length(t)) 'tps_' num2str(set) 'paramsetD'],'D')
                end
            end
        end
    end
end