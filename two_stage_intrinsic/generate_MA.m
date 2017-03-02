% This script generates the simulation files for the moment approximation.

clear all;
close all;
clc;

simMethod = 'MEC_2_LD_2_c';
%genSimFile('geneExp_MA','modelDef_toymodel',simMethod)

%%
k1 = 10; % basal mA/mB
k21 = 10; % mA stimulus
k22 = 20; % mB stimulus
k3 = 1; % deg mA/mB
k4 = 5; % A/B
k5 = 0.1; % deg A/B
t = linspace(0,60);

theta_s1 = [k1;k21;k3;k4;k5];
sol1 = simulate_geneExp_MA(t,theta_s1,[]);
theta_s2 = [k1;k22;k3;k4;k5];
sol2 = simulate_geneExp_MA(t,theta_s2,[]);

figure('name','MA')
l1 = plot(t,sol1.y(:,1),'-','color',color.SP); hold on;
fill([t fliplr(t)],...
    [sol1.y(:,1)'+sqrt(sol1.y(:,2))' fliplr(sol1.y(:,1)'-sqrt(sol1.y(:,2))')],1,...
    'facealpha',[0.1],'edgecolor','none', 'facecolor',color.SP); hold on;

l2 = plot(t,sol2.y(:,1),'--','color',color.SP); hold on;
fill([t fliplr(t)],...
    [sol2.y(:,1)'+sqrt(sol2.y(:,2))' fliplr(sol2.y(:,1)'-sqrt(sol2.y(:,2))')],1,...
    'facealpha',[0.1],'edgecolor','none', 'facecolor',color.SP); hold on;

xlabel('time [min]')
ylabel('A levels [au]')

legend([l1,l2],'subpopulation 1', 'subpopulation 2');
