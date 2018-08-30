clear all
close all
clc

load_plot_settings_robust

distributions = {'norm','skew_norm','students_t','neg_binomial'};
%%
options_plot.data.col{1} = color.data;
options_plot.data.bins = 40;
options_plot.model.lw =  1.5;
options_plot.boundaries(1).y_min = 0;
options_plot.boundaries(1).y_max = 250;
options_plot.subplot_lin = true; % subplots arranged linearly or as rectangle
options_plot.plainstyle =  true; % display axis
options_plot.legendflag =  false; % show legend
options_plot.titleflag = false; % show title
% Visualize representative data
for jDist = 1:length(distributions)
    for iDist = 1:length(distributions)
        if ~(iDist == 4 && jDist <4) % neg binomial only for neg binomial data
            load(['results/' distributions{jDist} '/results_geneExpression_simstudy_MA_' ...
                distributions{jDist} '_1000cells_3tps_1' distributions{iDist}],'M','D','parameters')
            
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
        end
    end
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.5 3])
    print('-depsc',['./figures/fits_geneExpression_MA_simstudy_' distributions{jDist}])
end

%% MSE and BIC

cellNumbers = [50,100,500,1000,5000];
timePoints{1} = [0,1,2];
timePoints{2} = [0,1,2,5];
timePoints{3} = [0,1,2,4,5];
timePoints{4} = [0,1,2,3,4,5];

confMatrix = zeros(4,4);
dist2true = nan(4,6,500);

count = 1;
for iDist = 1:4
    for iCells = 1:length(cellNumbers)
        for iTime = 1:length(timePoints)
            for iDataset = 1:5
                tmpBIC = inf(1,4);
                for jDist = 1:4
                    if ~(jDist == 4 && iDist <4) % neg binomial only for neg binomial data
                        clear parameters M D options conditions
                        load(['results/' distributions{iDist} '/results_geneExpression_simstudy_MA_' ...
                            distributions{iDist} ...
                            '_' num2str(cellNumbers(iCells)) 'cells_' ...
                            num2str(length(timePoints{iTime})) 'tps_' num2str(iDataset) ...
                            '' distributions{jDist}])
                        n_data = sum(sum(~isnan(D(1).y)));
                        parameters.MS.BIC = -2*parameters.MS.logPost(1) + ...
                            log(n_data)*parameters.number;
                        tmpBIC(jDist) = parameters.MS.BIC(1);
                    end
                    xi_true = D.xi_true;
                    if parameters.MS.par(6,1) > 0.5 % symmetry of mixture
                        parameters.MS.par(6,1) = 1-parameters.MS.par(6,1);
                        parameters.MS.par([2,3],1) = parameters.MS.par([3,2],1);
                    end
                    
                    dist2true(jDist,:,count) = parameters.MS.par(1:6,1)'-xi_true(1:6);
                end
                count = count+1;
                [~,indmin] = min(tmpBIC);
                confMatrix(iDist,indmin) = confMatrix(iDist,indmin)+1;
            end
        end
    end
end
%%
% remove neg binomial
confMatrix(1:3,4) = nan;
% Visualize confusion matrix
h=heatmap(confMatrix);
h.YDisplayLabels = {'normal','skew normal','Student''s t','negative binomial'};
h.XDisplayLabels = {'normal','skew normal','Student''s t','negative binomial'};
h.Colormap = flipud(color.greenwhite(30:end,2:4));
%flipud([zeros(64,1)        linspace(0,1,64)'  zeros(64,1);         % black to green
%linspace(0,1,64)'  ones(64,1)         linspace(0,1,64)']) ;%flipud(summer);
h.MissingDataColor = [0.8,0.8,0.8];
h.FontSize = 6;
h.ColorLimits = [0,100];
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7 5])
print('-depsc',['./figures/confusionmatrix'])

%%
figure('name','MSE SSA')
for iDist = 1:length(distributions)
    eval(['tmpcol = color.' distributions{iDist} ';']);
    meanMSE = nanmean(sum(dist2true(iDist,:,:).^2));
    stdMSE = nanstd(sum(dist2true(iDist,:,:).^2));
    plot(iDist, meanMSE, '.','MarkerSize',15,...
        'Color',tmpcol); hold on;
     plot([iDist,iDist], [meanMSE+stdMSE,meanMSE-stdMSE],'-','Color',tmpcol); hold on;
end
