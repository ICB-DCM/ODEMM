% This script generates the figures for the supplement for the conversion
% reaction example.

clear all
close all
clc

saveFigures = false;
load_plot_settings

%% plot all models
options_plot.data.col{1} = color.data;
options_plot.model.col{1} = color.SP;
models = {'RRE_onlyone','RRE_subpop','RRE_timedep','SP_k1', ...
    'SP_k1k2', 'SP_all', 'SP_k2', 'SP_k2k3', 'SP_k1k3', 'SP_k3'};

for m = 1:length(models)
    load(['./results/results_' models{m}],'parameters','M','D')
    if m > 3
        options_plot.simulate_musigma = 1; % for Sigma Points
    else
        options_plot.simulate_musigma = 0; % for RRE
    end
    xi = parameters.MS.par(:,1);
    fh = plotODEMM(D,M,xi,options_plot);
    set(fh{1},'name',['conversion reaction ' models{m}])
end

%% for reproducing figures of SI
options_plot.data.bins = 40;
options_plot.model.lw =  1;
options_plot.boundaries(1).x_min = 0.4; 
options_plot.boundaries(1).x_max = 0.8;
options_plot.boundaries(1).y_min = 0.35;
options_plot.boundaries(1).y_max = 0.83;
options_plot.subplot_lin = true; % subplots arranged linearly or as rectangle
options_plot.plainstyle =  true; % display axis
options_plot.legendflag =  false; % show legend
options_plot.titleflag = false; % show title

for m = 1:length(models)
    load(['./results/results_' models{m}],'parameters','M','D')
    if m > 3
        options_plot.simulate_musigma = 1; % for Sigma Points
    else
        options_plot.simulate_musigma = 0; % for RRE
    end
    xi = parameters.MS.par(:,1);
    plotODEMM(D,M,xi,options_plot);
    for i = 1:5
        s = subplot(1,5,i);
        ylim([0,0.15])
    end
    subplot(1,5,1);
    set(gca,'xtick',[0.4,0.8],...
        'xticklabel',{'0.4','0.8'},...
        'tickdir','out');
    ah=axes('position',[.13,.11,0.7813,0],'visible','off');
    line([0,1],[0,1],'parent',ah,'linewidth',0.6,'Color','k');
    ah=axes('position',[.133,.11,0,0.85],'visible','off');
    line([0,1],[0,1],'parent',ah,'linewidth',0.6,'Color','k');
    if saveFigures
        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4 2])
        set(gcf,'renderer','opengl');
        feval('print', '-depsc','-r1000',['./figures/datafit_cr_' num2str(m) '.eps']);
    end
end

%% Sampling vs Profiles
close all
load(['./results/results_SP_k3'],'parameters','M','D')
try
    load(['./results/results_BFchains_SP_k3'],'parameters_l')
    flagSampling = 1;
catch
    disp('did not find results_BFchains_SP_k3.mat. run main_sampling_ex1.m to obtain the samples.')
    flagSampling = 0;
end
thin = 1;
nbin = 10;
for iParam = 1:parameters.number
    subplot(3,3,iParam)
    if flagSampling
        [N,X] = hist(parameters_l{end}.S.par(iParam,parameters_l{end}.S.burnin+1:thin:end,1),nbin);
        h = bar(X,N/max(N),1,'facecolor',color.samp,'edgecolor',color.samp_edges); hold on;
    end
    plot(parameters.P(iParam).par(iParam,:),parameters.P(iParam).R,'LineWidth',1.2,...
        'color',color.opt);
    xlabel(parameters.name{iParam},'FontSize',6);
    box off
    set(gca,'TickDir','out')
end

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7 7.2])
    print('-depsc','./figures/sampleVSProfile');
end

%% Model selection
models = {'RRE_onlyone','RRE_subpop','RRE_timedep','SP_k1', 'SP_k1k2',...
    'SP_all', 'SP_k2', 'SP_k2k3', 'SP_k1k3', 'SP_k3'};

BICs = nan(10,1);
num_param = nan(10,1);
conv = nan(10,1);
lppds = nan(10,1);
logmargs = nan(10,1);
t_cpu_BF = zeros(10,2);
t_cpu_lppd = zeros(10,2);

for i = 1:10
    load(['./results/results_' models{i}],'parameters','D');
    BICs(i) = -2*parameters.MS.logPost(1)+ ...
        log(numel(D(1).t)*1000)*parameters.number;
    num_param(i) = parameters.number;
    t_cpu_opt(i) = nansum(parameters.MS.t_cpu(parameters.MS.t_cpu>0));
    try
        load(['./results/results_BFchains_' models{i}],'parameters_l');
        for iTemp = 1:17
            t_cpu_BF(i,1) = t_cpu_BF(i,1)+ parameters_l{iTemp}.S.t_cpu;
        end        
    catch
        disp('Sampling results not found.')
        t_cpu_BF(i,1) = NaN;
    end
    try
        clear parameters_l
        load(['./results/results_BFchains2_' models{i}],'parameters_l');
        for iTemp = 1:17
            t_cpu_BF(i,2) = t_cpu_BF(i,2)+ parameters_l{iTemp}.S.t_cpu;
        end
    catch
        disp('Sampling results not found.')
        t_cpu_BF(i,2) = NaN;
    end
    try
        load(['./results/lppd_results/results_lppd1_' models{i}],'parameters_lppd');
        for iData = 1:5
           t_cpu_lppd(i,1) = t_cpu_lppd(iData,1)+ parameters_lppd{iData}.S.t_cpu;
        end
        clear parameters_lppd
        load(['./results/lppd_results/results_lppd2_' models{i}],'parameters_lppd');
        for iData = 1:5
           t_cpu_lppd(i,2) = t_cpu_lppd(iData,2)+ parameters_lppd{iData}.S.t_cpu;
        end
    catch
        disp('Log pointwise predictive density results not found.')
        t_cpu_lppd(i,1) = NaN;
    end
end

% Plotting
figure
subplot('Position',[0.2,0.55,0.7,0.4]);

% lppd
load ./results/results_lppds
templppds = -lppds+max(max(lppds))+1;
plot(1:10,mean(templppds),'d','MarkerSize',3.5,'col',color.samp,'LineWidth',1); hold on;

% log marginals
load ./results/results_logmargs
tempmargs = -logmargs+max(max(logmargs))+1;
plot(1:10,nanmean(tempmargs),'.','MarkerSize',12,'col',color.samp); hold on;
plot([0,11],(log(100)+1)*ones(2,1),':','col',color.samp); hold on;

% BIC
plot(1:10,0.5*(BICs-min(BICs))+1,'.','MarkerSize',8,'col',color.opt); hold on;
plot([0,11],[6,6],':','col',color.opt); hold on;

box off
set(gca,'xtick',[1:10],...
    'ytick',[1,10,100,1000],...
    'xticklabel','',...
    'yscale','log',...
    'xlim',[0.5,10.5],...
    'ylim',[0.9,1001],...
    'TickDir','out','FontSize',6);

% correlations for the ranking provided by different model selection
% criteria
disp(['correlation BIC, lppd: ' num2str(spear(mean(templppds)',...
    BICs-min(BICs)+1))])
disp(['correlation BIC, BF: ' num2str(spear(mean(tempmargs)',...
    BICs-min(BICs)+1))])
disp(['correlation BF,lppd: ' num2str(spear(mean(templppds)',...
    mean(tempmargs)'))])

%% CPU time
subplot('Position',[0.2,0.1,0.7,0.4])
plot(1:10,t_cpu_opt,'.','MarkerSize',8,'col',color.opt,'LineWidth',1.1); hold on;
plot(1:10,nanmean(t_cpu_BF'),'.','MarkerSize',12,'col',color.samp); hold on;
plot(1:10,nanmean(t_cpu_lppd'),'d','MarkerSize',3.5,'col',color.samp,'LineWidth',1); hold on;
set(gca,...
    'xtick',[1:10],...
    'xticklabel','',...
    'ytick',[1e2,1e3,1e4,1e5],...
    'yscale','log',...
    'xlim',[0.5,10.5],...
    'ylim',[1e2,2e5],...
    'TickDir','out','FontSize',fs)
box off
if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4.5 6.3])
    print('-depsc','./figures/SIsamplingoptim')
end
%% Evidence for variability
figure
subplot('Position',[0.2,0.05,0.7,0.08])
diffs_b = BICs-min(BICs(1));
ws = exp(-0.5*diffs_b([1,4:end]));
ws = ws/nansum(ws);

probs(1) = sum(ws([1,2,3,6]+1));
probs(2) = sum(ws([2,3,4,6]+1));
probs(3) = sum(ws([3,5,6,7]+1));
load ./results/results_logmargs
for iRep = 1:size(logmargs,1)
    inds = [1,4:10];
    for i = 1:numel(inds)
        wssamp(i) = 1./(nansum(exp(logmargs(iRep,inds)-logmargs(iRep,inds(i)))));
    end
    probssamp(iRep,1) = nansum(wssamp([1,2,3,6]+1));
    probssamp(iRep,2) = nansum(wssamp([2,3,4,6]+1));
    probssamp(iRep,3) = nansum(wssamp([3,5,6,7]+1));
end

b = bar([probs;nanmean(probssamp,1)]','BarWidth',1); hold on;
b(1).FaceColor = color.opt;
b(2).FaceColor = color.samp;
box off
set(gca,'ytick',[0,0.5,1],...
    'xticklabel','',...
    'xlim',[0.65,3.35],...
    'TickDir','out','FontSize',fs);

plot_cr_variabilityReduction;

if saveFigures
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 3.5 10.5])
    print('-depsc','./figures/SI_cr_evidence')
end

