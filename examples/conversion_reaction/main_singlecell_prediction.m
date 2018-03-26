% This script shows the analysis of the single-cell trajectory prediction.
% First the single-cell parameter based on one single-cell trajectory are
% sampled and the predicted results are visualized. Then the parameters 
% are optimized for 100 single-cell trajectories. This analysis corresponds
% to Loos et al., Cell Systems (2018), Figure 3E-G.

clear all
close all
clc

% load population-based optimized parameters and the model
load('./results/results_SP_k3','parameters','M')

%% Example for sampling based on one single-cell trajectory calibrated  based on the last time point

% load single-cell trajectory
load('./data/sc_trajectory')

load_plot_settings
k1_1 = log(10)*parameters.MS.par(1,1);
k1_2 = log(10)*parameters.MS.par(1,2);
k2 = log(10)*parameters.MS.par(3,1);
k3 = theta(3);
m_k3 = log(10)*parameters.MS.par(4,1);
sigma_k3 = 10.^parameters.MS.par(5,1);
sigma_noise = 10.^parameters.MS.par(6,1);
w = parameters.MS.par(7,1);
k1 = k1_1;
clear parameters
parameters.name = {'k_1','k_3'};
parameters.number = 2;
parameters.min = [-3,-3];
parameters.max = [3,3];

% perform sampling
options.pesto = PestoOptions();
options.pesto.MCMC.nIterations         = 1e4;
options.pesto.MCMC.sigma0              = 10*diag(ones(1,parameters.number));
options.pesto.MCMC.theta0 = [-1;m_k3];
options.pesto.MCMC.samplingAlgorithm   = 'PT';
options.pesto.MCMC.PT.nTemps           = 1;
options.pesto.mode = 'text';
burnin = 1e3;
parameters = getParameterSamples(parameters, ...
    @(xi) sclogLikelihood_cr(t,ybar',xi,k2,m_k3,sigma_k3,sigma_noise,w,k1_1,k1_2,'last'), ...
    options.pesto);

% prediction for each sample
tsim = linspace(t(1),t(end));

ys = nan(numel(tsim),numel(burnin:options.pesto.MCMC.nIterations));
scount = 1;
for s = burnin:options.pesto.MCMC.nIterations
   k1 = parameters.S.par(1,s,1);
   if k1 < w
        k1 = k1_1;
    else
        k1 = k1_2;
   end
   sol = simulate_CR_log(tsim,[k1,k2,parameters.S.par(2,s,1)],[]);
   ys(:,scount) = sol.y;
   scount = scount +1;
end
ys = exp(ys);

% simulate true trajectory
sol = simulate_CR_log(tsim,theta,[]);
ytrue = sol.y;
%% Plot sampling results
figure
fill([tsim fliplr(tsim)],...
    [prctile(ys',99.5) ...
    fliplr(prctile(ys',0.5))],1,...
    'facealpha',[1],...
    'edgecolor','none',...
    'facecolor',color.CI_full(1,:)); hold on;

fill([tsim fliplr(tsim)],...
    [prctile(ys',80) ...
    fliplr(prctile(ys',10))],1,...
    'facealpha',[1],...
    'edgecolor','none',...
    'facecolor',color.CI_full(2,:)); hold on;

fill([tsim fliplr(tsim)],...
    [prctile(ys',55) ...
    fliplr(prctile(ys',45))],1,...
    'facealpha',[1],...
    'edgecolor','none',...
    'facecolor',color.CI_full(4,:)); hold on;

plot(tsim,median(ys'),'Color',color.SP); hold on;
plot(tsim,exp(ytrue),'-.','LineWidth',1.1,'Color',color.true); hold on;
plot(tsim(end),ybar(end),'.','MarkerSize',10,'Color','k'); hold on;

% formatting figure
set(gca,'TickDir','out',...
    'xtick',[0,1,2],'xticklabel',[0,60,120],'FontSize',6)
set(gca,'TickDir','out',...
    'ytick',[0.5,0.6,0.7],'xticklabel',[0,60,120],'FontSize',6)
set(gca,'TickDir','out',...
    'ytick',[0.5,0.6,0.7],'xticklabel',[0,60,120],'FontSize',6,'ylim',[0.48,0.72])
box off
xlabel('time [min]','FontSize',6)
ylabel('B levels [au]','FontSize',6)

%% save figure
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4 3.2])
%set(gcf,'renderer','opengl');
%feval('print', '-depsc','-r1000',['./figures/SI_cr_scpred_lastttp.eps']);
%%
i = 1;
param{i} = parameters;
param{i}.guess = [k1_1,k1_2;m_k3,m_k3];
param{i} =  getMultiStarts(param{i}, ...
        @(xi) sclogLikelihood_cr(t,ybar',xi,k2,m_k3,sigma_k3,sigma_noise,w,k1_1,k1_2,'last'), ...
        options.MS);

%% Optimization for all single-cell trajectories
load sc_trajectories_100
options.MS = PestoOptions();
options.MS.mode = 'silent';
options.MS.n_starts = 20;
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions  = optimset('algorithm','interior-point',...
    'MaxIter',5000,...
    'display','off',...
    'GradObj','on',...
    'Hessian','off',...
    'TolFun',1e-10, 'TolX',1e-10,...
    'PrecondBandWidth',inf,...
    'Display','iter-detailed',...
    'MaxFunEvals',1e5);

for i = 1:100
    param{i} = parameters;
    param{i}.guess = [k1_1,k1_2;m_k3,m_k3];
    param{i} =  getMultiStarts(param{i}, ...
        @(xi) sclogLikelihood_cr(t,ybar(:,i)',xi,k2,m_k3,sigma_k3,sigma_noise,w,k1_1,k1_2,'last'), ...
        options.MS);
    
   k1 = param{i}.MS.par(1,1);
   if k1 < w
        k1 = k1_1;
    else
        k1 = k1_2;
   end
   sol = simulate_CR_log([0 0.1 0.5 1 2],[k1,k2,param{i}.MS.par(2,1)],[]);
   ypred(:,i) = exp(sol.y);
end
%%
close all
figure
subplot(1,3,1)
hs=scatter(ypred(1,:),y(1,:),30,'.'); hold on;
set(hs,'MarkerEdgeColor',color.Erk);
set(hs,'MarkerEdgeAlpha',0.5);
title(['r=' num2str(corr(ypred(1,:)',y(1,:)'),3)],'FontSize',6)
axis square
ylabel('true B levels','FontSize',6)
xlabel('predicted B levels','FontSize',6)
box off
xlim([0.4,0.8]); ylim([0.4,0.8]);
set(gca,'TickDir','out')

subplot(1,3,2)
hs2=scatter(ypred(3,:),y(3,:),30,'.'); hold on;
set(hs2,'MarkerEdgeColor',color.k5);
set(hs2,'MarkerEdgeAlpha',0.5);   
title(['r=' num2str(corr(ypred(3,:)',y(3,:)'),3)],'FontSize',6)
axis square
xlabel('predicted B levels','FontSize',6)
box off
xlim([0.4,0.8]); ylim([0.4,0.8]);
set(gca,'TickDir','out')

subplot(1,3,3)
hs2=scatter(ypred(end,:),y(end,:),30,'.'); hold on;
set(hs2,'MarkerEdgeColor',color.k1);
set(hs2,'MarkerEdgeAlpha',0.5);
title(['r=' num2str(corr(ypred(end,:)',y(end,:)'),3)],'FontSize',6)
axis square
xlabel('predicted B levels','FontSize',6)
box off
set(gca,'TickDir','out')
xlim([0.4,0.8]); ylim([0.4,0.8]);

%% save figure
%set(gca,'TickDir','Out')
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7.5 2.3])
%print('-dpdf','./figures/SI_crprediction_scatter');
