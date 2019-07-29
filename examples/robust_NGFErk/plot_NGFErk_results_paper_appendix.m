% Plot model fits for full data set
close all
clear all
clc

%% Plot only 2D-data (kernel density estimate)
distribution = 'students_t';
load_plot_settings_robust
load(['./results/results_NGFErk_' distribution],'parameters');
load('./data/log_data_PDL');
n_data = 0;
for e = 1:4
    n_data = n_data + sum(sum(sum(~isnan(D(e).y))));
end

eval(['ODEMM_NGFErk_' distribution]);

options_plot.x_scale = 'lin';
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.data.lw = 1;
options_plot.model.lw = 1.2;
options_plot.model.level_linewidth = 0.8;
options_plot.data.col = color.data;
options_plot.boundaries(1).y_min = 4; options_plot.boundaries(1).y_max = 10;
options_plot.boundaries(2).y_min = 4; options_plot.boundaries(2).y_max = 10;
options_plot.boundaries(3).y_min = [1*log(10),2*log(10)]; % 1 x n_meas vector
options_plot.boundaries(3).y_max = log([2e4,2e4]); % 1 x n_meas vector
options_plot.boundaries(4).y_min = log([9e1,1e2]); % 1 x n_meas vector
options_plot.boundaries(4).y_max = log([5e3,2e4]);

for d = 1:6
    options_plot.model.levelsets{4,d} = [0.15 0.35 0.65 0.9];
    options_plot.model.levelsets{3,d} = [0.08 0.16 0.32 0.48];
end
options_plot.subplot_lin = true;
options_plot.marginals = false;
options_plot.data.markersize = 10;
options_plot.data.kde = true;
options_plot.plainstyle = true;
options_plot.titleflag = false;
options_plot.legendflag = false;
options_plot.I = [3,4];
options_plot.tu_ind{3} = [1:6];
options_plot.tu_ind{4} = [1:6];
options_plot.model.col = 'k';
fh = plotODEMM(D,[],[],options_plot);
for e = 3:4
    tmpf = figure(fh{e});
    subplot(1,length(options_plot.tu_ind{e}),1)
    set(gca,'xtick',[2*log(10),3*log(10),4*log(10)],...
        'xticklabel',{'10^2','10^3','10^4'});
    set(gca,'ytick',[2*log(10),3*log(10),4*log(10)],...
        'yticklabel',{'10^2','10^3','10^4'});
    set(gca,'FontSize',8)
    for s = 2:length(options_plot.tu_ind{e})
        subplot(1,length(options_plot.tu_ind{e}),s)
        set(gca,'yticklabel','')
        set(gca,'xtick',[2*log(10),3*log(10),4*log(10)],...
        'xticklabel',{'10^2','10^3','10^4'});
        set(gca,'FontSize',8)
    end
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 2.5])
    print('-dpng',['figures/data_NGFErk_e' num2str(e) '_all'],'-r1000')
end

%% Plot model fits for 2D-data
distributions = {'norm','skew_norm','students_t'};
for iDist = 1:numel(distributions)
    options_plot.I = [3,4];
    load(['./results/results_NGFErk_' distributions{iDist}],'parameters');
    eval(['ODEMM_NGFErk_' distributions{iDist}]);
    eval(['options_plot.model.col = color.' distributions{iDist}]);
    options_plot.data.kde = false;
    options_plot.marginals = true;
    xi = parameters.MS.par(:,1);
    [conditions,D] = collectConditions(D,M);
    fh = plotODEMM(D,M,xi,options_plot);
    for e = 3:4
        tmpf = figure(fh{e});
        subplot(1,length(options_plot.tu_ind{e}),1)
        set(gca,'ytick',[2*log(10),3*log(10),4*log(10)],'yticklabel',{'10^2','10^3','10^4'});
        set(gca,'xtick',[2*log(10),3*log(10),4*log(10)],'xticklabel',{'10^2','10^3','10^4'});
        set(gca,'FontSize',8)
        for s = 2:length(options_plot.tu_ind{e})
            subplot(1,length(options_plot.tu_ind{e}),s)
            set(gca,'ytick',[2*log(10),3*log(10),4*log(10)],'yticklabel',{''});
            set(gca,'xtick',[2*log(10),3*log(10),4*log(10)],'xticklabel',...
                {'10^2','10^3','10^4'});
            set(gca,'FontSize',8)
        end
        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 2.5])
        print('-dpng',['figures/fit_NGFErk ' distributions{iDist} '_e'...
            num2str(e) '_all'],'-r1000')
    end   
end

%% Plot model fits for 1D-data
options_plot.tu_ind{1} = [1:7];
for iDist = 1:numel(distributions)
    load(['./results/results_NGFErk_' distributions{iDist}],'parameters');
    eval(['ODEMM_NGFErk_' distributions{iDist}]);
    eval(['options_plot.model.col = color.' distributions{iDist}]);
    eval(['options_plot.model.ls = linestyles.' distributions{iDist}]);
    eval(['options_plot.model.lw = linewidths.' distributions{iDist} ';']);
    options_plot.data.kde = false;
    xi = parameters.MS.par(:,1);
    [conditions,D] = collectConditions(D,M);
    options_plot.I = [1];
    options_plot1 = options_plot;
    Dtmp = D;
    if iDist > 1
        options_plot1.fh = fh1;
        options_plot1.hold_on = true;
        Dtmp(1).y(:) = nan;
    end
    fh1 = plotODEMM(Dtmp,M,xi,options_plot1);  
    
end
set(gca,'xtick',[2*log(10),3*log(10),4*log(10)],'xticklabel',...
    {'10^2','10^3','10^4'},'TickDir','out');
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 6.5 2.5])
print('-dpdf',['figures/fit_NGFErk_e1_all'])

options_plot.tu_ind{2} = [1:7];
for iDist = 1:numel(distributions)
    load(['./results/results_NGFErk_' distributions{iDist}],'parameters');
    eval(['ODEMM_NGFErk_' distributions{iDist}]);
    eval(['options_plot.model.col = color.' distributions{iDist}]);
    eval(['options_plot.model.ls = linestyles.' distributions{iDist}]);
    eval(['options_plot.model.lw = linewidths.' distributions{iDist} ';']);
    options_plot.data.kde = false;
    xi = parameters.MS.par(:,1);
    [conditions,D] = collectConditions(D,M);
    options_plot.I = [2];
    options_plot2 = options_plot;
    Dtmp = D;
    if iDist > 1
        options_plot2.fh = fh2;
        options_plot2.hold_on = true;
        Dtmp(2).y(:) = nan;
    end
    fh2 = plotODEMM(Dtmp,M,xi,options_plot2); 
    tmpf = figure(fh2{2});
end
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 6.5 2.5])
print('-dpdf',['figures/fit_NGFErk_e2_all'])

