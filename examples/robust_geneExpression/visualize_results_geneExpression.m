clear all
close all
clc

load_plot_settings_robust

%% SSA results
load data/data_geneExpression_SSA
n_data = sum(sum(~isnan(D(1).y)));
xi_true = D(1).xi_true;

distributions = {'norm','skew_norm','students_t','neg_binomial'};
distribution_names = {'normal','skew normal','Student''s t','negative binomial'};

BICs = nan(numel(distributions),1);
llhs = nan(numel(distributions),1);
dist2trues = nan(numel(distributions),numel(xi_true));

for iDist = 1:length(distributions)
    load(['./results/results_geneExpression_SSA_' distributions{iDist}],'parameters')
    BICs(iDist) = -2*parameters.MS.logPost(1) + log(n_data)*parameters.number;
    llhs(iDist) = parameters.MS.logPost(1);
    parameters.MS.par(numel(xi_true),1)
    if parameters.MS.par(numel(xi_true),1) > 0.5 % symmetrie of mixture
        parameters.MS.par(numel(xi_true),1) = 1-parameters.MS.par(numel(xi_true),1);
        parameters.MS.par([2,3],1) = parameters.MS.par([3,2],1);
    end
    dist2trues(iDist,:) = parameters.MS.par(1:numel(xi_true),1)'-xi_true;
end

figure('name','geneExpression SSA')
for iDist = 1:length(distributions)
    subplot(2,1,1);
    eval(['tmpcol = color.' distributions{iDist} ';']);
    plot(iDist, BICs(iDist)-min(BICs), '.','MarkerSize',15,'Color',tmpcol); hold on;
    subplot(2,1,2);
    plot(iDist, llhs(iDist), '.','MarkerSize',15,'Color',tmpcol); hold on;
end
subplot(2,1,1);
box off
set(gca,'xtick',[1:length(distributions)],'xticklabels','',...
    'TickDir','out')
ylabel('\Delta BIC')
subplot(2,1,2);
box off
set(gca,'xtick',[1:length(distributions)],'xticklabels','',...
    'TickDir','out')
ylabel('log-likelihood')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7 5])
print('-depsc',['./figures/BICllh_geneExpression_SSA'])

%%
figure('name','MSE SSA')
for iDist = 1:length(distributions)
    eval(['tmpcol = color.' distributions{iDist} ';']);
    plot(iDist, sum(dist2trues(iDist,:).^2), '.','MarkerSize',15,...
        'Color',tmpcol); hold on;
end
box off
set(gca,'xtick',[1:length(distributions)],'xticklabels','',...
    'TickDir','out')
ylabel('MSE')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7 2])
print('-depsc',['./figures/MSE_geneExpression_SSA'])

%% Plot fitting results
options_plot.data.col{1} = color.data;
options_plot.data.bins = 40;
options_plot.model.lw =  1.5;
options_plot.boundaries(1).y_min = 0; 
options_plot.boundaries(1).y_max = 2000;
options_plot.subplot_lin = true; % subplots arranged linearly or as rectangle
options_plot.plainstyle =  true; % display axis
options_plot.legendflag =  false; % show legend
options_plot.titleflag = false; % show title
for iDist = 1:length(distributions)
    load(['./results/results_geneExpression_SSA_' distributions{iDist}])
    eval(['tmpcol = color.' distributions{iDist} ';']);
    eval(['tmpls = linestyles.' distributions{iDist} ';']);
    eval(['tmplw = linewidths.' distributions{iDist} ';']);

    options_plot.model.col{1} = tmpcol;
    options_plot.model.ls =  tmpls;
    options_plot.model.lw =  tmplw;
    xi = parameters.MS.par(:,1);
    if iDist == 1
        options_plot.hold_on = false;
        fh = plotODEMM(D,M,xi,options_plot);
    else
        options_plot.hold_on = true;
        options_plot.fh = fh;
        D(1).y(:) = nan;
        fh = plotODEMM(D,M,xi,options_plot);
    end
    for iSubplot = 1:length(D.t)
        subplot(1,length(D.t),iSubplot)
        if(iSubplot==1)
            set(gca,'FontSize',6,'TickDir','out')
        end
        ylim([0,0.35]);
    end
%     set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 10 4])
%     print('-depsc',['./figures/fit_geneExpression_SSA_' distributions{iDist}])
%     %print('-dpdf',['./figures/fit_geneExpression_SSA_' distributions{iDist} '.pdf'])
end
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 11 4])
print('-depsc',['./figures/fits_geneExpression_SSA'])

%% 
for jDist = 1%:length(distributions)
    for iDist = 1:length(distributions)
        load(['./results/results_geneExpression_MA_' distributions{jDist} ...
            '_' distributions{iDist}])
        eval(['tmpcol = color.' distributions{iDist} ';']);
        eval(['tmpls = linestyles.' distributions{iDist} ';']);
        eval(['tmplw = linewidths.' distributions{iDist} ';']);

        options_plot.model.col{1} = tmpcol;
        options_plot.model.ls =  tmpls;
        options_plot.model.lw =  tmplw;
        xi = parameters.MS.par(:,1);
        if iDist == 1
            options_plot.hold_on = false;
            fh = plotODEMM(D,M,xi,options_plot);
        else
            options_plot.hold_on = true;
            options_plot.fh = fh;
            D(1).y(:) = nan;
            fh = plotODEMM(D,M,xi,options_plot);
        end
        for iSubplot = 1:length(D.t)
            subplot(1,length(D.t),iSubplot)
            if(iSubplot==1)
                set(gca,'FontSize',6,'TickDir','out')
            end
            ylim([0,0.4]);
        end
    end
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 10 3])
    print('-depsc',['./figures/fits_geneExpression_MA_' distributions{jDist}])
end
