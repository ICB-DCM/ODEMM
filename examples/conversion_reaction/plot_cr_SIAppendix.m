% This script generates the figures for the supplement

clear all
close all
clc
load_plot_settings
%% Figure SI Appendix
options_plot.x_scale = 'lin';
options_plot.data.plot = 'filled';
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.boundaries(1).x_min = [0.4];
options_plot.boundaries(1).x_max = [0.8];
options_plot.subplot_lin = 1;
options_plot.plainstyle =  0;
options_plot.legendflag =  1;
options_plot.model.lw =  1;
options_plot.data.lw =  0.7;
options_plot.data.fill_col{1} = color.data;
options_plot.model.col{1} = color.SP;
options_plot.boundaries(1).y_min = 0.35;
options_plot.boundaries(1).y_max = 0.83;

models = {'RRE_onlyone','RRE_subpop','RRE_timedep','SP_k1', 'SP_k1k2', 'SP_all', 'SP_k2', 'SP_k2k3', 'SP_k1k3', 'SP_k3'};
for m = 1:10
    load(['./project/results/results_' models{m}],'parameters','M','D')
    if m > 3
        options_plot.simulate_musigma = 1; % for Sigma Points
    else
        options_plot.simulate_musigma = 0; % for RRE
    end
    xi = parameters.MS.par(:,1);
    plotODEMix(D,M,xi,[],options_plot);
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
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4 2])
    set(gcf,'renderer','opengl');
    feval('print', '-depsc','-r1000',['./project/figures/datafit_cr_' num2str(m) '.eps']);
end

%% Sampling vs Profiles
close all
load(['./project/results/results_SP_k3_profiles'],'parameters','M','D')
load(['./project/results/results_BayesFactors/results_BFchains_SP_k3'],'parameters_l')
thin = 1;
nbin = 10;
for iParam = 1:parameters.number
    %     if iParam == parameters.number
    %         subplot(3,3,9)
    %     else
    subplot(3,3,iParam)
    % end
    [N,X] = hist(parameters_l{end}.S.par(iParam,parameters_l{end}.S.burnin+1:thin:end,1),nbin);
    h = bar(X,N/max(N),1,'facecolor',color.samp,'edgecolor',color.samp_edges); hold on;
    plot(parameters.P(iParam).par(iParam,:),parameters.P(iParam).R,'LineWidth',1.2,...
        'color',color.opt);
    xlabel(parameters.name{iParam},'FontSize',6);
    box off
    set(gca,'TickDir','out')
end
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7 7.2])
print('-depsc','./project/figures/sampleVSProfile');


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
    load(['./project/results/results_' models{i}],'parameters','D');
    BICs(i) = -2*parameters.MS.logPost(1)+ ...
        log(numel(D(1).t)*1000)*parameters.number;
    num_param(i) = parameters.number;
    t_cpu_opt(i) = nansum(parameters.MS.t_cpu(parameters.MS.t_cpu>0));
    try
        load(['./project/results/results_BayesFactors/results_BFchains_' models{i}],'parameters_l');
        t_cpu_samp(i,1) = parameters_l{end}.S.t_cpu;
        for iTemp = 1:17
            t_cpu_BF(i,1) = t_cpu_BF(i,1)+ parameters_l{iTemp}.S.t_cpu;
        end        
    catch
         %t_cpu_samp(i,1) = NaN;
         t_cpu_BF(i,1) = NaN;
    end
    try
        clear parameters_l
        load(['./project/results/results_BayesFactors/results_BFchains2_' models{i}],'parameters_l');
        %t_cpu_samp(i,2) = parameters_l{end}.S.t_cpu;
        for iTemp = 1:17
            t_cpu_BF(i,2) = t_cpu_BF(i,2)+ parameters_l{iTemp}.S.t_cpu;
        end
    catch
         t_cpu_samp(i,2) = NaN;
         t_cpu_BF(i,2) = NaN;
    end
    try
        load(['./project/results/lppd_results/results_lppd1_' models{i}],'parameters_lppd');
        for iData = 1:5
           t_cpu_lppd(i,1) = t_cpu_lppd(iData,1)+ parameters_lppd{iData}.S.t_cpu;
        end
        clear parameters_lppd
        load(['./project/results/lppd_results/results_lppd2_' models{i}],'parameters_lppd');
        for iData = 1:5
           t_cpu_lppd(i,2) = t_cpu_lppd(iData,2)+ parameters_lppd{iData}.S.t_cpu;
        end
    catch
    end
end
%% BIC
figure
subplot('Position',[0.2,0.55,0.7,0.4]);
% lppd
load ./project/results/lppd_results/results_lppds
templppds = -lppds+max(max(lppds))+1;
plot(1:10,mean(templppds),'d','MarkerSize',3.5,'col',color.samp,'LineWidth',1); hold on;
% for m = 1:10
%     myerrorbar(m,mean(templppds(:,m)'),std(templppds(:,m)),0.07,color.samp_edges);
% end
% log marginal likelihood
load ./project/results/results_BayesFactors/results_logmargs
logmargs = logmargs';
tempmargs = -logmargs+max(max(logmargs))+1;
plot(1:10,nanmean(tempmargs),'.','MarkerSize',12,'col',color.samp); hold on;
% for m = 1:10
%     myerrorbar(m,mean(tempmargs(:,m)'),std(tempmargs(:,m)),0.07,color.samp_edges);
% end
plot([0,11],(log(100)+1)*ones(2,1),':','col',color.samp); hold on;
% BIC
plot(1:10,0.5*(BICs-min(BICs))+1,'.','MarkerSize',8,'col',color.opt); hold on;
plot([0,11],[6,6],':','col',color.opt); hold on;

box off
set(gca,'xtick',[1:10],...
    'ytick',[1,10,100,1000],...'yminortick','off'
    ...'yticklabel',{''10^{1}'','10^{1}','10^2','10^2','10^3'},...
    'xticklabel','',...
    'yscale','log',...
    'xlim',[0.5,10.5],...
    'ylim',[0.9,1001],...
    'TickDir','out','FontSize',6);

disp(['correlation BIC, lppd: ' num2str(spear(mean(templppds)',...
    BICs-min(BICs)+1))])
disp(['correlation BIC, BF: ' num2str(spear(mean(tempmargs)',...
    BICs-min(BICs)+1))])
disp(['correlation BF,lppd: ' num2str(spear(mean(templppds)',...
    mean(tempmargs)'))])
%% CPU
subplot('Position',[0.2,0.1,0.7,0.4])
plot(1:10,t_cpu_opt,'.','MarkerSize',8,'col',color.opt,'LineWidth',1.1); hold on;

%plot(1:10,mean(t_cpu_samp'),'+','MarkerSize',3,'col',color.samp,'LineWidth',1.1); hold on;
% for m = 1:10
%     myerrorbar(m,mean(t_cpu_samp(m,:)'),std(t_cpu_samp(m,:)'),0.07,color.samp);
% end

plot(1:10,nanmean(t_cpu_BF'),'.','MarkerSize',12,'col',color.samp); hold on;
% for m = 1:10
%     myerrorbar(m,mean(t_cpu_BF(m,:)'),std(t_cpu_BF(m,:)'),0.07,color.samp_edges);
% end

plot(1:10,mean(t_cpu_lppd'),'d','MarkerSize',3.5,'col',color.samp,'LineWidth',1); hold on;
% for m = 1:10
%     myerrorbar(m,mean(t_cpu_lppd(m,:)'),std(t_cpu_lppd(m,:)'),0.07,color.samp_edges);
% end
set(gca,...
    ...'ytick',[1,11,101,1001],...
    ...'yticklabel',{'0','10^1','10^2','10^3'},...
    'xtick',[1:10],...
    'xticklabel','',...
    'ytick',[1e2,1e3,1e4,1e5],...
    'yscale','log',...
    'xlim',[0.5,10.5],...
    'ylim',[1e2,2e5],...
    'TickDir','out','FontSize',fs)
box off
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4.5 6.3])
print('-depsc','./project/figures/SIsamplingoptim')

%% Evidence for variability
figure
subplot('Position',[0.2,0.05,0.7,0.08])
diffs_b = BICs-min(BICs(1));
ws = exp(-0.5*diffs_b([1,4:end]));
ws = ws/nansum(ws);

probs(1) = sum(ws([1,2,3,6]+1));
probs(2) = sum(ws([2,3,4,6]+1));
probs(3) = sum(ws([3,5,6,7]+1));
load ./project/results/results_BayesFactors/results_logmargs
logmargs = logmargs';

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
% myerrorbar(1.15,mean(probssamp(1,:)'),std(probssamp(1,:)),0.07);
% myerrorbar(2.15,mean(probssamp(2,:)'),std(probssamp(2,:)),0.07);
% myerrorbar(3.15,mean(probssamp(3,:)'),std(probssamp(3,:)),0.07);

box off
set(gca,'ytick',[0,0.5,1],...
    'xticklabel','',...
    'xlim',[0.65,3.35],...
    'TickDir','out','FontSize',fs);
%%
plot_variabilityReduction;
%%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 3.5 10.5])
print('-depsc','./project/figures/SI_crevidence')

