% Visualization script for the fits

close all
clear all
clc

load ./project/data/data_Lys_1D2D_sepscaled
%load ./project/results/results_NGF_subpop_final_T
load_plot_settings
%%
%ODEMM_NGF_subpop_T_final; % TrkA

TextSizes.DefaultAxesFontSize = 6;
TextSizes.DefaultTextFontSize = 6;
set(0,TextSizes);

options_plot.x_scale = 'log';
options_plot.data.plot = 'filled';
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.data.lw = 1;
options_plot.model.lw = 0.8;
options_plot.data.fill_col{1} = color.Lys_data;
options_plot.data.fill_col{2} = color.Lys_data;
options_plot.data.fill_col{3} = color.Lys_data;
options_plot.data.fill_col{4} = color.Lys_data;

options_plot.data.col{1} = color.Lys_data;
options_plot.model.col{1} =  color.Lys;
options_plot.data.col{2} = options_plot.data.col{1};
options_plot.model.col{2} =  options_plot.model.col{1};
options_plot.data.col{3} = options_plot.data.col{2};
options_plot.data.col{4} = options_plot.data.col{2};
options_plot.model.col{3} = options_plot.model.col{1};
options_plot.model.col{4} = options_plot.model.col{1};


for d = 1:6
    options_plot.model.levelsets{4,d} = [0.9,1.7,2.5,3.3];
    options_plot.model.levelsets{3,d} = [0.45 0.75 1.15 1.5];
end
options_plot.model.levelsets{4,1} = [1.3,2.6,4,5];
options_plot.model.levelsets{4,3} = [0.9,1.7,2.5,4];
options_plot.fit_scaling_parameters = true;
options_plot.boundaries(1).y_min = 1e2; options_plot.boundaries(1).y_max = 2e4;
options_plot.boundaries(2).y_min = 1e2; options_plot.boundaries(2).y_max = 2e4;
options_plot.boundaries(3).y_min = [1e1,1e2]; % 1 x n_meas vector
options_plot.boundaries(3).y_max = [2e4,2e4]; % 1 x n_meas vector
options_plot.boundaries(4).y_min = [9e1,1e2]; % 1 x n_meas vector
options_plot.boundaries(4).y_max = [5e3,2e4]; % 1 x n_meas vector
options_plot.subplot_lin = 1;

options_plot.switch_axes = true;
options_plot.marginals = 0;
options_plot.sameplot = 0;
options_plot.model.subpopulations = 0;
options_plot.simulate_musigma = true;
options_plot.data.markersize = 10;
options_plot.model.level_linewidth = 0.5;
tu_ind{3} = [1:6];
tu_ind{4} = [1:6];
%xi = parameters.MS.par(:,1);
%[conditions,D] = collectConditions(D,M);
close all
%[fh] = plotODEMix(D,M,xi,[3:4],options_plot,tu_ind);
[fh] = plotODEMix(D,[],[],[3:4],options_plot,tu_ind);

mtvals(1:8) = log10([2:9]*10^1);
mtvals(9:16) = log10([2:9]*10^2);
mtvals(17:24) = log10([2:9]*10^3);
%% TrkA/pErk
figure(fh{3});

for i = 6:-1:1
    s=subplot(1,6,i)
    view(180,90);
    set(gca,'xdir','reverse');
    set(gca,'ytick', [2,3,4]);
    
    hA = gca;
    hA.YAxis.MinorTickValues = mtvals;
    hA.XAxis.MinorTickValues = mtvals;
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(s, 'position', [0.1+(i-1)*0.14,0.2,0.13,0.25] );
    
    set(gca, 'xtick', [2,3,4] ,'yticklabel',{''},...
        'xticklabel',{'10^2','','10^4'});
    
    
    if i == 1
        set(gca  ,'yticklabel',{'10^2','10^3','10^4'});
    end
end
%%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 9 4.5])
print('-dpdf','./project/figures/data_TrkApErk_all')
%% Erk/pErk
figure(fh{4});
for i = 6:-1:1
    s=subplot(1,6,i)
    view(180,90);
    set(gca,'xdir','reverse');
    set(gca, 'xtick', [2,3] ,'xticklabel',{''},'ytick', [2,3,4],...
        'yticklabel',{''});
    hA = gca;
    hA.YAxis.MinorTickValues = mtvals;
    hA.XAxis.MinorTickValues = mtvals;
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(s, 'position', [0.1+(i-1)*0.14,0.5,0.13,0.25] );
    
    set(gca, 'xtick', [2,3] ,'yticklabel',{''},...
        'xticklabel',{'10^2','10^3'});
    
    
    if i == 1
        set(gca  ,'yticklabel',{'10^2','10^3','10^4'});
    end
end
%%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 9 4.5])
print('-dpdf','./project/figures/data_ErkpErk_all')


close all

