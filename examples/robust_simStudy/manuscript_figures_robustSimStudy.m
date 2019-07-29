% Generates the boxplots and example data histogram.

clear all
close all
clc

noiseFlag = true;
modelnames = {'conversionReaction','diffProteinExpression',...
    'twoStageGeneExpression','diffProteinExpression'};
distributions = {'norm','skew_norm','students_t','neg_binomial'};
outlierstrs = {'outlier2_zeros','outlier5_dublets','outlier10_unif'};
load_simStudy_settings

for m = [1,3,4]
    BICs = [];
    llhs = [];
    MSEs = [];
    converged = [];
    t_cpus = [];
    all_t_cpus{1} = [];
    all_t_cpus{2} = [];
    all_t_cpus{3} = [];
    all_t_cpus{4} = [];
    
    %% No outliers
    count = 1;
    for it = 1:3
        t = tps{it};
        for ic = 1:4
            for set = 1:3
                for iDist = 1:4
                    try
                        if noiseFlag && m == 3 && iDist == 4 && ...
                                strcmp(outlierstr,'outlier5_dublets')
                            load (['./resultsNoise/results_noise_model' ...
                                num2str(m) '_twoStage_' num2str(n_cells(ic)) ...
                                'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                                'paramset_' distributions{iDist}],'D','parameters');
                        elseif noiseFlag
                            load (['./resultsNoise/results_noise_model' ...
                                num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                                'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                                'paramset_' distributions{iDist}],'D','parameters');
                        else
                            load (['./results/results_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                                'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                                'paramset_' distributions{iDist}],'D','parameters');
                        end
                        switch m
                            case 1
                                ind_weight = 5;
                                ind_s1 = 2;
                                ind_s2 = 1;
                            case 2
                                ind_weight = 7;
                                ind_s1 = [3,5];
                                ind_s2 = [2,4];
                            case 3
                                ind_weight = 7;
                                ind_s1 = 3;
                                ind_s2 = 2;
                            case 4
                                ind_weight = 5;
                                ind_s1 = [2];
                                ind_s2 = [3];
                        end
                        if D.xi_true(ind_weight) == 0.5 % check which MSE is lower)
                            D.xi_true2 = D.xi_true;
                        end
                        % wrongly assigned in generation of data
                        D.xi_true([ind_s1,ind_s2]) = D.xi_true([ind_s2,ind_s1]);
                        % symmetry of mixture
                        if ~(sign(parameters.MS.par(ind_weight,1)-0.5) == sign(D.xi_true(ind_weight)-0.5))
                            parameters.MS.par(ind_weight,1) = 1-parameters.MS.par(ind_weight,1);
                            parameters.MS.par([ind_s1,ind_s2],1) = parameters.MS.par([ind_s2,ind_s1],1);
                        end
                        allmodelnames{count} = ['model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                            'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                            'paramset'];
                        MSEs(count,iDist) = sum((D.xi_true'-parameters.MS.par(1:length(D.xi_true),1)).^2);
                        if D.xi_true(ind_weight) == 0.5 % check which MSE is lower)
                            MSE2 = sum((D.xi_true2'-parameters.MS.par(1:length(D.xi_true),1)).^2);
                            MSEs(count,iDist) = min(MSEs(count,iDist),MSE2);
                        end
                        BICs(count,iDist) = -2*parameters.MS.logPost(1) + ...
                            log(sum(sum(sum(~isnan(D(1).y)))))*(parameters.number);
                        llhs(count,iDist) = parameters.MS.logPost(1);
                        converged(count,iDist) = sum(parameters.MS.logPost(1) - parameters.MS.logPost<1e-3);
                        t_cpus(count,iDist) = sum(parameters.MS.t_cpu);
                        all_t_cpus{iDist} = [all_t_cpus{iDist};parameters.MS.t_cpu];
                    catch e
                        disp(e.message)
                        MSEs(count,iDist) = NaN;
                        BICs(count,iDist) = NaN;
                        llhs(count,iDist) = NaN;
                        converged(count,iDist) = NaN;
                        t_cpus(count,iDist) = NaN;
                        all_t_cpus{iDist} = [all_t_cpus{iDist}; nan(30,1)];
                    end
                end
                count = count+1;
            end
        end
    end
    fh = robustSimStudy_boxplots(BICs,MSEs,all_t_cpus,converged,t_cpus);
    figure(fh)
    fh.PaperPosition=[0 0 4 9];
    if noiseFlag
        print('-depsc',['./figures/nooutlier_boxplot_' modelnames{m}]);
    else
        print('-depsc',['./figures/nooutlier_noise_boxplot_' modelnames{m}]);
    end
    %% All outlier scenarios merged
    BICs = [];
    llhs = [];
    MSEs = [];
    converged = [];
    t_cpus = [];
    all_t_cpus{1} = [];
    all_t_cpus{2} = [];
    all_t_cpus{3} = [];
    all_t_cpus{4} = [];
    count = 1;
    for iOut = 1:3
        outlierstr = outlierstrs{iOut};
        for it = 1:3
            t = tps{it};
            for ic = 1:4
                for set = 1:3
                    for iDist = 1:length(distributions)
                        try
                            if noiseFlag
                                load(['./resultsOutlierNoise/results_noise_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                                    'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                                    'paramset_' distributions{iDist} '_' outlierstr],'D','parameters');
                            else
                                load(['./resultsOutlier/results_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                                    'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                                    'paramset_' distributions{iDist} '_' outlierstr],'D','parameters');
                            end
                            switch m
                                case 1
                                    ind_weight = 5;
                                    ind_s1 = 2;
                                    ind_s2 = 1;
                                case 2
                                    ind_weight = 7;
                                    ind_s1 = [3,5];
                                    ind_s2 = [2,4];
                                    ind_s2 = 2;
                                case 4
                                    ind_weight = 5;
                                    ind_s1 = [2];
                                    ind_s2 = [3];
                            end
                            if D.xi_true(ind_weight) == 0.5 % check which MSE is lower)
                                D.xi_true2 = D.xi_true;
                            end
                            % wrongly assigned in generation of data
                            D.xi_true([ind_s1,ind_s2]) = D.xi_true([ind_s2,ind_s1]);
                            % symmetry of mixture
                            if ~(sign(parameters.MS.par(ind_weight,1)-0.5) == sign(D.xi_true(ind_weight)-0.5))
                                parameters.MS.par(ind_weight,1) = 1-parameters.MS.par(ind_weight,1);
                                parameters.MS.par([ind_s1,ind_s2],1) = parameters.MS.par([ind_s2,ind_s1],1);
                            end
                            allmodelnames{count} = ['model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                                'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                                'paramset_' outlierstr];
                            MSEs(count,iDist) = sum((D.xi_true'-parameters.MS.par(1:length(D.xi_true),1)).^2);
                            BICs(count,iDist) = -2*parameters.MS.logPost(1) + ...
                                log(sum(sum(sum(~isnan(D(1).y)))))*(parameters.number);
                            llhs(count,iDist) = parameters.MS.logPost(1);
                            converged(count,iDist) = sum(parameters.MS.logPost(1) - parameters.MS.logPost<1e-3);
                            t_cpus(count,iDist) = sum(parameters.MS.t_cpu);
                            all_t_cpus{iDist} = [all_t_cpus{iDist};parameters.MS.t_cpu];
                        catch e
                            disp(['/results_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
                                'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                                'paramset_' distributions{iDist} '_' outlierstr])
                            disp(e.message)
                            MSEs(count,iDist) = NaN;
                            BICs(count,iDist) = NaN;
                            llhs(count,iDist) = NaN;
                            converged(count,iDist) = NaN;
                            t_cpus(count,iDist) = NaN;
                            all_t_cpus{iDist} = [all_t_cpus{iDist}; nan(30,1)];
                        end
                    end
                    count = count+1;
                end
            end
        end
    end
    % Generate boxplots
    fh = robustSimStudy_boxplots(BICs,MSEs,all_t_cpus,converged,t_cpus);
    figure(fh)
    fh.PaperPosition=[0 0 4 9];
    if noiseFlag
        print('-depsc',['./figures/outlier_noise_boxplot_' modelnames{m}]);
    else
        print('-depsc',['./figures/outlier_boxplot_' modelnames{m}]);
    end
end

%% Outlier histograms
close all
clear options_plot
load_plot_settings_robust
m = 1;
ic = 4;
it = 3; t = tps{it};
set = 1;
positions{1} = [0.1,0.2,0.25,0.7];
positions{2} = [0.4,0.2,0.25,0.7];
positions{3} = [0.7,0.2,0.25,0.7];
load(['./data/data_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
    'cells_' num2str(length(t)) 'tps_' num2str(set) 'paramsetD'],'D')
y_orig = squeeze(D.y(1,end,:));
outlierstrs = {'outlier2_zeros','outlier5_dublets','outlier10_unif'};
for iOut = 1:3
    load(['./dataOutlier/data_model' num2str(m) '_' modelnames{m} '_' num2str(n_cells(ic)) ...
        'cells_' num2str(length(t)) 'tps_' num2str(set) 'paramset_' outlierstrs{iOut}],'D')
    D5 = D;
    ind_nooutlier5 = find(y_orig==squeeze(D5.y(1,end,:)));
    ind_outlier5 = find(~(y_orig==squeeze(D5.y(1,end,:))));
    switch iOut
        case 1
            histEdges = linspace(0,400,51)';
        case 2
            histEdges = linspace(150,800,51)';
        case 3
            histEdges = linspace(150,450,51)';
    end
    y_out5 = squeeze(D5.y(1,end,ind_outlier5));
    y_noout5 =  squeeze(D5.y(1,end,ind_nooutlier5));
    fhtmp{iOut}=figure;
    ax = axes('Parent',fhtmp{iOut});
    hold(ax,'on');
    ylim([0,300]);
    ax.FontSize = fs;
    ax.TickDir='out';
    
    h = hist(y_noout5,histEdges);
    h=h(1:end-1)';
    noOut.data = fill(histEdges(round(0.5:0.5:length(histEdges))),...
        [0;h(round(0.5:0.5:length(h)));0], color.data); hold on;
    noOut.data.EdgeColor = color.data;
    noOut.data.FaceAlpha = 1;
    hold on;
    
    hOut = hist(y_out5,histEdges);
    hOut=hOut(1:end-1)';
    noOut.data = fill(histEdges(round(0.5:0.5:length(histEdges))),...
        [0;hOut(round(0.5:0.5:length(hOut)));0], color.outlier); hold on;
    noOut.data.EdgeColor = color.outlier;
    noOut.data.FaceAlpha = 1;
    plot([0,1000],[0,0],'k-','LineWidth',0.5)
    box off
    switch iOut
        case 1
            xlim([0,400]);
        case 2
            xlim([100,800]);
        case 3
            xlim([150,450]);
    end
    xlabel('')
    ylabel('')
    fhtmp{iOut}.PaperPosition=[0 0 6 4];
    print('-depsc',['./figures/' outlierstrs{iOut} '_histogram']);
end

