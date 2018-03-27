% This script visualizes the results for the conversion reaction.

clear all
close all
clc

load_plot_settings

%% Figure 2C
load('./results/results_SP_k3','parameters','M','D')
xi = parameters.MS.par(:,1);

options_plot.simulate_musigma = 1; % required for simulation
options_plot.data.col{1} = color.data;
options_plot.model.col{1} = color.SP;
options_plot.data.bins = 40;
plotODEMM(D,M,xi,options_plot);

%% to reproduce figure for paper
options_plot.model.lw =  1;
options_plot.boundaries(1).x_min = 0.4; 
options_plot.boundaries(1).x_max = 0.8;
options_plot.boundaries(1).y_min = 0.35;
options_plot.boundaries(1).y_max = 0.83;
options_plot.subplot_lin = true; % subplots arranged linearly or as rectangle
options_plot.plainstyle =  true; % display axis
options_plot.legendflag =  false; % show legend
options_plot.titleflag = false; % show title
plotODEMM(D,M,xi,options_plot);
for i = 1:5
    s = subplot(1,5,i);
    ylim([0,0.15])
end
subplot(1,5,1);
set(gca,'xtick',[0.4,0.8],'xticklabel',{'0.4','0.8'});
set(gca,'tickdir','out')
ah=axes('position',[.13,.11,0.7813,0],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.6,'Color','k');
ah=axes('position',[.133,.11,0,0.85],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.6,'Color','k');
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 10 4])
%print('-dpdf','./figures/data_cr.pdf')

%% Figure 2B - BIC
models = {'RRE_onlyone','RRE_subpop','RRE_timedep','SP_k1', 'SP_k1k2',...
    'SP_all', 'SP_k2', 'SP_k2k3', 'SP_k1k3', 'SP_k3'};

t_cpus = nan(10,100);
BICs = nan(10,1);
num_param = nan(10,1);
conv = nan(10,1);
for i = 1:10
    load(['./results/results_' models{i}],'parameters','D');
    BICs(i) = -2*parameters.MS.logPost(1)+ log(numel(D(1).t)*1000)*parameters.number;
    num_param(i) = parameters.number;
    t_cpus(1,1:numel(parameters.MS.t_cpu(parameters.MS.t_cpu>0))) = parameters.MS.t_cpu(parameters.MS.t_cpu>0);
    conv(i) = sum(2*(parameters.MS.logPost-parameters.MS.logPost(1))>-icdf('chi2',0.95,1))/parameters.MS.n_starts;
end

% BIC
subplot('Position',[0.2,0.3,0.7,0.25]);
plot(1:10,BICs-min(BICs)+1,'.','MarkerSize',8,'col','k'); hold on;
ylabel('BIC-min(BIC)','FontSize',6);
set(gca,'yminortick','off',...
    'ytick',[1,11,101,1001],...
    'yticklabel',{'0','10^1','10^2','10^3'},...
    'xtick',[1:10],...
    'xticklabel','',...
    'yscale','log',...
    'xlim',[0.5,10.5],...
    'ylim',[0.7,1001],...
    'TickDir','out')
box off

% parameter number
subplot('Position',[0.2,0.6,0.7,0.25]);
plot(1:10,num_param,'.','MarkerSize',8,'col','k'); hold on;
ylabel({'number of','parameters'},'FontSize',6)
set(gca,'xtick',[1:10],...
    'xticklabel','',...
    'xlim',[0.5,10.5],...
    'ylim',[4,16],...
    'TickDir','out');
box off
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4.5 5])
%print('-depsc','./figures/BICs')

%% Confidence intervals variability
alpha = [0.99,0.95,0.9,0.8];
delta = 0;
epsilon1 = 2e-3;
epsilon2 = 2e-3;
xi_true = [-0.1,0.1,-0.45,-0.2,-1,-1.8,0.7];

load('./results/results_SP_k3')
options.MS.mode = 'silent';
figure('name','confidence intervals');
lh = [];
xi_min = -2;
xi_max = 0.6;

subplot(2,1,1)
for j = 1:length(alpha)
    pars = getParameterConfidenceIntervals(parameters,alpha(j),options.MS);
    pars.CI.PL = max(min(pars.CI.PL,xi_max),xi_min);
    fill(pars.CI.PL(5,[1,2,2,1]),((-2-delta))+(0.95+j*0.05)*0.5*[+1,+1,-1,-1],...
        'b','facecolor',color.SP_CI(j,:)); hold on;
    plot([xi_true(5), xi_true(5)],[-3,-1],'-','color',color.true,'LineWidth',0.8);
end
plot([pars.MS.par(5),pars.MS.par(5)],((-2-delta))+0.7*[+1,-1],...
    'LineWidth',1,'color',color.SP_CI(end,:)); hold on;
box off
set(gca,'xscale','lin',...
    'xlim',[-1.05,-0.95],...
    'xtick',[-1.05,-1,-0.95],'xticklabel',{'10^{-1.05}','10^{-1}', '10^{-0.95}'},...
    'ytick',[],...
    'Tickdir','out',...
    'yticklabel',{});

subplot(2,1,2)
for j = 1:length(alpha)
    pars = getParameterConfidenceIntervals(parameters,alpha(j),options.MS);
    pars.CI.PL = max(min(pars.CI.PL,xi_max),xi_min);
    fill(pars.CI.PL(6,[1,2,2,1]),((-2-delta))+(0.95+j*0.05)*0.5*[+1,+1,-1,-1],...
        'b','facecolor',color.SP_CI(j,:)); hold on;
    plot([xi_true(6), xi_true(6)],[-3,-1],'-','color',color.true,'LineWidth',0.8);    
end
plot([pars.MS.par(6),pars.MS.par(6)],((-2-delta))+0.7*[+1,-1],...
    'LineWidth',1,'color',color.SP_CI(end,:)); hold on;
box off
set(gca,'xscale','lin',...
    'xlim',[-1.9,-1.5],...
    'xtick',[-1.9,-1.7,-1.5],'xticklabel',{'10^{-1.9}','10^{-1.7}', '10^{-1.5}'},...
    'ytick',[],...
    'Tickdir','out',...
    'yticklabel',{});
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 3 2])
%print('-dpdf','./figures/uncertainty_mainman');

% add labels 
subplot(2,1,1)
ylabel('\sigma_{k_3}')
subplot(2,1,2)
ylabel('\sigma_{noise}')
xlabel('parameter value')
