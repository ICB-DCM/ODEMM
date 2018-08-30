clear all
close all
clc

load_plot_settings_robust

%%%%%%%%%%%%%%%%%%%% SSA %%%%%%%%%%%%%%%%%%%
%% BIC llh and MSE
load data/data_geneExpression_2_SSA
n_data = sum(sum(~isnan(D(1).y)));
xi_true = D(1).xi_true;

distributions = {'norm','skew_norm','students_t','neg_binomial'};
distribution_names = {'normal','skew normal','Student''s t','negative binomial'};

BICs = nan(numel(distributions),1);
llhs = nan(numel(distributions),1);
dist2trues = nan(numel(distributions),numel(xi_true));

for iDist = 1:length(distributions)
    load(['./results/results_geneExpression_2_SSA_' distributions{iDist}],'parameters')
    xi_true
    BICs(iDist) = -2*parameters.MS.logPost(1) + log(n_data)*parameters.number;
    llhs(iDist) = parameters.MS.logPost(1);
    parameters.MS.par(numel(xi_true),1)
    if parameters.MS.par(numel(xi_true),1) > 0.5 % symmetry of mixture
        parameters.MS.par(numel(xi_true),1) = 1-parameters.MS.par(numel(xi_true),1);
        parameters.MS.par([2,3],1) = parameters.MS.par([3,2],1);
    end
    dist2trues(iDist,:) = parameters.MS.par(1:numel(xi_true),1)'-xi_true;
end

%%
figure('name','geneExpression SSA')
for iDist = 1:length(distributions)
    subplot(2,1,1);
    eval(['tmpcol = color.' distributions{iDist} ';']);
    plot(iDist, BICs(iDist)-min(BICs), '.','MarkerSize',15,'Color',tmpcol); hold on;
    subplot(2,1,2);
    plot(iDist, -llhs(iDist), '.','MarkerSize',15,'Color',tmpcol); hold on;
end
subplot(2,1,1);
box off
set(gca,'xtick',[1:length(distributions)],'xticklabels','',...
    'TickDir','out')
ylabel('\Delta BIC','FontSize',bigfs)
subplot(2,1,2);
box off
set(gca,'xtick',[1:length(distributions)],'xticklabels','',...
    'TickDir','out')
ylabel('neg. log-likelihood','FontSize',bigfs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7 5.5])
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
ylabel('MSE','FontSize',bigfs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7 2])
print('-depsc',['./figures/MSE_geneExpression_SSA'])

%%
figure('name','geneExpression SSA all')
for iDist = 1:length(distributions)
    subplot(3,1,1);
    eval(['tmpcol = color.' distributions{iDist} ';']);
    plot(iDist, BICs(iDist)-min(BICs), '.','MarkerSize',15,'Color',tmpcol); hold on;
    subplot(3,1,2);
    plot(iDist, -llhs(iDist), '.','MarkerSize',15,'Color',tmpcol); hold on;
    subplot(3,1,3);
    plot(iDist, sum(dist2trues(iDist,:).^2), '.','MarkerSize',15,...
        'Color',tmpcol); hold on;
end
subplot(3,1,1);
box off
set(gca,'xtick',[1:length(distributions)],'xticklabels','',...
    'TickDir','out')
ylabel('\Delta BIC','FontSize',bigfs)
subplot(3,1,2);
box off
set(gca,'xtick',[1:length(distributions)],'xticklabels','',...
    'TickDir','out')
ylabel({'negative'; 'log-likelihood'},'FontSize',bigfs)
subplot(3,1,3);
box off
set(gca,'xtick',[1:length(distributions)],'xticklabels','',...
    'TickDir','out')
ylabel('MSE','FontSize',bigfs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7 6])
print('-depsc',['./figures/BICllhMSE_geneExpression_SSA'])

%% Profiles 
load data/data_geneExpression_2_SSA
n_data = sum(sum(~isnan(D(1).y)));
xi_true = D(1).xi_true;
distributions = {'norm','skew_norm','students_t','neg_binomial'};

arrange2x3 = false;

for iParam = 1:numel(xi_true)
    if arrange2x3
    	subplot(2,3,iParam);
    else
        subplot(1,6,iParam);
    end
    for iDist = 1:length(distributions)
        load(['./results/results_geneExpression_2_SSA_' distributions{iDist}],'parameters')
        eval(['tmpcol = color.' distributions{iDist} ';']);
        eval(['tmpls = linestyles.' distributions{iDist} ';']);
        eval(['tmplw = linewidths.' distributions{iDist} ';']);
        symmetryChange = false;
        if parameters.MS.par(numel(xi_true),1) > 0.5 % symmetry of mixture
            symmetryChange = true;
        end
        
        if iDist == 2
            plot([xi_true(iParam), xi_true(iParam)],[0,1],'Color',color.true,...
                'LineWidth',2); hold on;
        end
        
        if iParam == 6 && symmetryChange
            plot(1-parameters.P(6).par(6,:),parameters.P(6).R,...
                'Color',tmpcol,'LineStyle',tmpls,'LineWidth',tmplw);
        elseif iParam == 2 && symmetryChange
            plot(parameters.P(3).par(3,:),parameters.P(3).R,...
                'Color',tmpcol,'LineStyle',tmpls,'LineWidth',tmplw);
        elseif iParam == 3 && symmetryChange
            plot(parameters.P(2).par(2,:),parameters.P(2).R,...
                'Color',tmpcol,'LineStyle',tmpls,'LineWidth',tmplw);
        else
            plot(parameters.P(iParam).par(iParam,:),parameters.P(iParam).R,...
                'Color',tmpcol,'LineStyle',tmpls,'LineWidth',tmplw);
        end
    end
    xlabel(parameters.name{iParam},'FontSize',bigfs)
    box off
    set(gca,'TickDir','out')
    if ~ismember(iParam,[1,4]) && arrange2x3
        %set(gca,'yticklabel','');
    elseif iParam>1 && ~arrange2x3
    	%set(gca,'yticklabel','');
    else
        ylabel('likelihood ratio','FontSize',bigfs)
    end
end
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 20 3])
print('-depsc',['./figures/profiles_geneExpression_SSA'])

%% Plot fitting results
options_plot.data.col{1} = color.data;
options_plot.data.bins = 40;
options_plot.model.lw =  1.5;
options_plot.boundaries(1).y_min = 0;
options_plot.boundaries(1).y_max = 300;
options_plot.subplot_lin = true; % subplots arranged linearly or as rectangle
options_plot.plainstyle =  true; % display axis
options_plot.legendflag =  false; % show legend
options_plot.titleflag = false; % show title
for iDist = 1:length(distributions)
    load(['./results/results_geneExpression_2_SSA_' distributions{iDist}])
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
        ylim([0,0.25]);
    end
    %     set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 10 4])
    %     print('-depsc',['./figures/fit_geneExpression_SSA_' distributions{iDist}])
    %     %print('-dpdf',['./figures/fit_geneExpression_SSA_' distributions{iDist} '.pdf'])
end
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 15 4])
print('-depsc',['./figures/fits_geneExpression_SSA'])

%%
%%%%%%%%%%%%%%%%% Moment Approximation results %%%%%%%%%%%%%%%%%%%%%
for jDist = 1:length(distributions)
    for iDist = 1:length(distributions)
        load(['./results/results_geneExpression_2_MA_' distributions{jDist} ...
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
