% This script visualizes the results for the conversion reaction

clear all
close all
clc

load_plot_settings
%% Figure 2C
load('./project/results/results_SP_k3','parameters','M','D')
options_plot.simulate_musigma = 1;

xi = parameters.MS.par(:,1);
options_plot.x_scale = 'lin';
options_plot.data.plot = 'filled';
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.boundaries(1).x_min = [0.4]; % 1 x n_meas vector
options_plot.boundaries(1).x_max = [0.8];
options_plot.subplot_lin = 1;
options_plot.model.lw =  1;
options_plot.data.lw =  0.7;
options_plot.data.fill_col{1} = color.data;
options_plot.model.col{1} = color.SP;
options_plot.model.ls = '-';

options_plot.boundaries(1).y_min = 0.35;
options_plot.boundaries(1).y_max = 0.83;
%plotODEMix(D,M,xi,[],options_plot);
plotODEMix(D,[],[],[],options_plot);
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
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4 2])
set(gcf,'renderer','opengl');
feval('print', '-depsc','-r1000','./project/figures/data_cr.eps');
%%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 10 4])
print('-dpdf','./project/figures/data_cr.pdf')
%% Figure 2B - BIC

models = {'RRE_onlyone','RRE_subpop','RRE_timedep','SP_k1', 'SP_k1k2', 'SP_all', 'SP_k2', 'SP_k2k3', 'SP_k1k3', 'SP_k3'};

t_cpus = nan(10,100);
BICs = nan(10,1);
num_param = nan(10,1);
conv = nan(10,1);
for i = 1:10
    load(['./project/results/results_' models{i}],'parameters','D');
    BICs(i) = -2*parameters.MS.logPost(1)+ log(numel(D(1).t)*1000)*parameters.number;
    num_param(i) = parameters.number;
    t_cpus(1,1:numel(parameters.MS.t_cpu(parameters.MS.t_cpu>0))) = parameters.MS.t_cpu(parameters.MS.t_cpu>0);
    conv(i) = sum(2*(parameters.MS.logPost-parameters.MS.logPost(1))>-icdf('chi2',0.95,1))/parameters.MS.n_starts;
end

% BIC
%BICs = BICs(sort_ind);
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
print('-depsc','./project/figures/BICs')

%% Confidence intervals variability
alpha = [0.99,0.95,0.9,0.8];
delta = 0;
epsilon1 = 2e-3;
epsilon2 = 2e-3;
xi_true = [-0.1,0.1,-0.45,-0.2,-1,-1.8,0.7];

load('./project/results/results_SP_k3')
figure('name','Confidence intervals');
lh = [];
xi_min = -2;
xi_max = 0.6;

subplot(2,1,1)
for j = 1:length(alpha)
    pars = getParameterConfidenceIntervals(parameters,alpha(j));
    pars.CI.PL = max(min(pars.CI.PL,xi_max),xi_min);
    fill(pars.CI.PL(5,[1,2,2,1]),((-2-delta))+(0.95+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.SP_CI(j,:)); hold on;
    %lh(1) = fill((pars.MS.par(5)+epsilon1*0.15*[-1,+1,+1,-1]),((-2-delta))+0.7*[+1,+1,-1,-1],'b','facecolor',color.SP_CI(end,:)); hold on;
    plot([xi_true(5), xi_true(5)],[-3,-1],'-','color',color.true,'LineWidth',0.8);
end
plot([pars.MS.par(5),pars.MS.par(5)],((-2-delta))+0.7*[+1,-1],'LineWidth',1,'color',color.SP_CI(end,:)); hold on;
box off
set(gca,'xscale','lin',...
    'xlim',[-1.05,-0.95],...
    'xtick',[-1.05,-1,-0.95],'xticklabel',{'10^{-1.05}','10^{-1}', '10^{-0.95}'},......
    ...'xtick',[log10(0.095),log10(0.1),log10(0.105)],...
    ...'xticklabel',{'0.095','0.1', '0.015'},...
    'ytick',[],...
    'Tickdir','out',...
    'yticklabel',{});

subplot(2,1,2)
for j = 1:length(alpha)
    pars = getParameterConfidenceIntervals(parameters,alpha(j));
    pars.CI.PL = max(min(pars.CI.PL,xi_max),xi_min);
    fill(pars.CI.PL(6,[1,2,2,1]),((-2-delta))+(0.95+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.SP_CI(j,:)); hold on;
    %lh(1) = fill((pars.MS.par(6)+epsilon1*1*[-1,+1,+1,-1]),((-2-delta))+0.7*[+1,+1,-1,-1],'b','facecolor',color.SP_CI(end,:)); hold on;
    plot([xi_true(6), xi_true(6)],[-3,-1],'-','color',color.true,'LineWidth',0.8);    
end
plot([pars.MS.par(6),pars.MS.par(6)],((-2-delta))+0.7*[+1,-1],'LineWidth',1,'color',color.SP_CI(end,:)); hold on;
box off
set(gca,'xscale','lin',...
    'xlim',[-1.9,-1.5],...
    'xtick',[-1.9,-1.7,-1.5],'xticklabel',{'10^{-1.9}','10^{-1.7}', '10^{-1.5}'},......
    ...'xtick',log10([0.012 0.02 0.028]),...
    ...'xticklabel',{'0.012','0.02', '0.028'},...
    'ytick',[],...
    'Tickdir','out',...
    'yticklabel',{});
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 3 2])
print('-dpdf','./project/figures/uncertainty_mainman');

%% Supplement - Convergence
% figure
% conv = conv([1,2,3,   4     7    10     5     9     8     6]);
%
% b = bar([conv;ones(1,numel(conv))]);
%
% for i = 1:numel(conv)
%     if i < 4
%         b(i).FaceColor = color.RRE;
%     else
%         b(i).FaceColor = color.SP;
%     end
% end
% xlim([0.6,1.4])
% set(gca,'xtick',[]);
% set(gca,'ytick',[0,0.5,1],'yticklabel',{'0','50','100'});
% set(gca,'TickDir','out');
% ylabel('convergence [%]')
% box off;
%
% set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4.5 3])
% print('-dpdf','./project/figures/convergence_RREvsSP');

%% Supplement - tcpu
% b = boxplot(t_cpus');
% % for i = 1:3
% % set(boxpl(i), 'Color', col_RRE);
% % end
% % for i = 4:10
% % set(boxpl(i), 'Color', col_SP);
% % end
% set(gca,'TickDir','out');
% set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4 3])
% print('-dpdf','./project/figures/tcpu_RREvsSP');