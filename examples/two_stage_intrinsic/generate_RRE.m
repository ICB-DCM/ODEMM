% This script generates the simulation files for the reaction rate
% equations
%
% @note The simulation file for the RRE needs to be modified manually after
% generation, use the already generated and modified file
% RRE_geneExp_RRE_syms. The simulation file from that needs to be generated
% using cvodwrap.

clear all;
close all;
clc;

%% file generation
% simMethod = 'RRE';
% genSimFile('geneExp_RRE','modelDef_toymodel',simMethod)
cvodewrap('geneExp_RRE','RRE_geneExp_RRE_syms');
%% 
k1 = 10; % basal mA/mB
k21 = 10; % mA stimulus
k22 = 20; % mB stimulus
k3 = 1; % deg mA/mB
k4 = 5; % A/B
k5 = 0.1; % deg A/B
t = linspace(0,60);

theta_s1 = [k1;k21;k3;k4;k5];
sol1 = simulate_geneExp_RRE(t,theta_s1,[]);
theta_s2 = [k1;k22;k3;k4;k5];
sol2 = simulate_geneExp_RRE(t,theta_s2,[]);

figure('name','RRE')
plot(t,sol1.y(:,1),'-','color',color.RRE); hold on;
plot(t,sol2.y(:,1),'--','color',color.RRE); hold on;
xlabel('time [min]')
ylabel('A levels [au]')
legend('subpopulation 1', 'subpopulation 2');
