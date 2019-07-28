clear all
close all
clc

outlierstrs = {'outlier2_zeros','outlier5_dublets','outlier10_unif'};
noiseFlag = true;

for iOut = 1:3
    close all
    clear options_plot
    outlierstr = outlierstrs{iOut};
    tps{1}=[0,0.5,1,2,4];
    tps{2}=[0,0.5,2,4];
    tps{3}=[0,0.5,2];
    n_cells = [50,100,500,1000];

    modelnames = {'conversionReaction','diffProteinExpression','twoStageGeneExpression'};
    distributions = {'norm','neg_binomial','skew_norm','students_t'};

    load_plot_settings_robust
    m = 1;
    ic = 4;
    it = 3; t = tps{it};
    set = 1;

    warning off
    load(['./dataOutlier/data_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
        'cells_' num2str(length(t)) 'tps_' num2str(set) 'paramset_' outlierstr],'D')
    warning on
    options_plot.x_scale = 'lin';
    options_plot.data.bins = 50;
    options_plot.model.points = 200;
    options_plot.data.lw = 1;
    options_plot.model.lw = 1.2;
    options_plot.model.level_linewidth = 0.8;
    options_plot.data.col = color.data;
    options_plot.subplot_lin = false;
    options_plot.marginals = false;
    options_plot.data.markersize = 10;
    options_plot.data.kde = true;
    options_plot.plainstyle = false;
    options_plot.titleflag = false;
    options_plot.legendflag = false;
    options_plot.y_counts = true;
    % plot only data
    options_plot.I = [1];
    options_plot.tu_ind{1} = [3];
    switch outlierstr
        case 'outlier2_zeros'
            options_plot.boundaries(1).y_min = 0;
            options_plot.boundaries(1).y_max = 400;
        case 'outlier10_unif'
            options_plot.boundaries(1).y_min = 150;
            options_plot.boundaries(1).y_max = 450;
        case 'outlier5_dublets'
            options_plot.boundaries(1).y_min = 150;
            options_plot.boundaries(1).y_max = 800;
    end
    options_plot.model.col = 'k';
    Dtmp = D;
    for iDist = 1:numel(distributions)
        if noiseFlag
            load(['./resultsOutlierNoise/results_noise_model' num2str(m) '_' ...
                modelnames{m} '_' num2str(n_cells(ic)) ...
                 'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                 'paramset_' distributions{iDist} '_' outlierstr]);
        else
        load(['./resultsOutlier/results_model' num2str(m) '_' modelnames{m}...
            '_' num2str(n_cells(ic)) ...
            'cells_' num2str(length(t)) 'tps_' num2str(set) ...
            'paramset_' distributions{iDist} '_' outlierstr]);
        end
        eval(['options_plot.model.col = color.' distributions{iDist}]);
        eval(['options_plot.model.ls = linestyles.' distributions{iDist}]);
        eval(['options_plot.model.lw = linewidths.' distributions{iDist}]);
        xi = parameters.MS.par(:,1);
        [conditions,D] = collectConditions(D,M);
        if iDist > 1
            options_plot.fh = fh;
            options_plot.hold_on = true;
            Dtmp(1).y(:) = nan;
        end
        fh = plotODEMM(Dtmp,M,xi,options_plot);
    end
    box off
    ax = get(gca);
    ylim([0,0.3]);
    ax.TickDir = 'out';
    ax.FontSize = fs;
    xlabel('')
    ylabel('')
    fh{1}.PaperPosition=[0 0 6 4];
    if noiseFlag
    	print('-depsc',['./figures/' outlierstr '_fits_noise']);
    else
        print('-depsc',['./figures/' outlierstr '_fits']);
    end
end
