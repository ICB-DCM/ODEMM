close all

load('./project/data/data_matrices_1D2D');

load ./project/results/results_4_1802

p = parameters{15};
clear parameters;
parameters = p;
ODEMM_NGF_ECM_T_comb15

%% Visualization - Fit
load_plot_settings
options_plot.x_scale = 'log';
options_plot.data.plot = 'filled';

options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.data.lw = 1;
options_plot.model.lw = 1;

options_plot.data.col{1} = color.Lys_data;
options_plot.model.col{1} =  color.Lys;
options_plot.model.col{3} =  color.Lys;
options_plot.model.col{5} =  color.Lys;
options_plot.model.col{7} =  color.Lys;

options_plot.data.col{2} = color.ColI_data;
options_plot.model.col{2} = color.ColI;
options_plot.model.col{4} =  color.ColI;
options_plot.model.col{6} =  color.ColI;
options_plot.model.col{8} =  color.ColI;
options_plot.data.col{3} = options_plot.data.col{1};
options_plot.data.col{5} = options_plot.data.col{1};
options_plot.data.col{7} = options_plot.data.col{1};
options_plot.data.col{4} = options_plot.data.col{2};
options_plot.data.col{6} = options_plot.data.col{2};
options_plot.data.col{8} = options_plot.data.col{2};
% options_plot.model.colormap{6} =  'autumn';
% options_plot.model.colormap{8} =  'autumn';
% options_plot.model.colormap{5} =  'winter';
% options_plot.model.colormap{7} =  'winter';

options_plot.model.levelsets{8,1} = [1.3,2.6,4,5];
options_plot.model.levelsets{8,3} = [0.9,1.7,2.5,4];
options_plot.model.levelsets{8,6} = [0.9,1.7,2.5,3.3];
options_plot.model.levelsets{7,1} = [1.3,2.6,4,5];
options_plot.model.levelsets{7,3} = [0.9,1.7,2.5,4];
options_plot.model.levelsets{7,6} = [0.9,1.7,2.5,3.3];

options_plot.model.levelsets{5,1} = [0.45 0.75 1.15 1.5];
options_plot.model.levelsets{5,3} = [0.45 0.75 1.15 1.5];
options_plot.model.levelsets{5,6} = [0.45 0.75 1.15 1.5];
options_plot.model.levelsets{6,1} = [0.45 0.75 1.15 1.5];
options_plot.model.levelsets{6,3} = [0.45 0.75 1.15 1.5];
options_plot.model.levelsets{6,6} = [0.45 0.75 1.15 1.5];

options_plot.data.fill_col{3} = options_plot.data.col{3};
options_plot.data.fill_col{4} = options_plot.data.col{4};
options_plot.data.fill_col{5} = options_plot.data.col{5};
options_plot.data.fill_col{6} = options_plot.data.col{6};
options_plot.data.fill_col{7} = options_plot.data.col{7};
options_plot.data.fill_col{8} = options_plot.data.col{8};

options_plot.fit_scaling_parameters = true;
options_plot.boundaries(1).y_min = 1e2; options_plot.boundaries(1).y_max = 2e4;
options_plot.boundaries(2).y_min = 1e2; options_plot.boundaries(2).y_max = 2e4;
options_plot.boundaries(3).y_min = 1e2; options_plot.boundaries(3).y_max = 2e4;
options_plot.boundaries(4).y_min = 1e2; options_plot.boundaries(4).y_max = 2e4;

options_plot.boundaries(5).y_min = [1e0,5e1];
options_plot.boundaries(5).y_max = [1e4,5e5];
options_plot.boundaries(6).y_min = [1e0,5e1];
options_plot.boundaries(6).y_max = [1e4,5e5];
options_plot.boundaries(7).y_min = [1e2,1e2];
options_plot.boundaries(7).y_max = [5e3,2e5];
options_plot.boundaries(8).y_min = [1e2,1e2];
options_plot.boundaries(8).y_max = [5e3,2e5];

options_plot.switch_axes = true;
options_plot.marginals = 1;
options_plot.sameplot = 1;
options_plot.model.subpopulations = 0;
options_plot.simulate_musigma = true;
D(1).measurand = 'pErk conc.';
D(2).measurand = 'pErk conc.';
D(3).measurand = 'pErk conc.';
D(4).measurand = 'pErk conc.';
%tu_ind{1} = [1:7];
%tu_ind{2} = [1:7];
tu_ind{3} = [1,3,7];
tu_ind{4} = [1,3,7];
tu_ind{5} = [1];
tu_ind{6} = [1];
tu_ind{7} = [1];
tu_ind{8} = [1];
xi = parameters.MS.par(:,1);
[conditions,D] = collectConditions(D,M);
%%
close all
fh = plotODEMix_ECM(D,M,xi,[3,4,5,6,7,8],options_plot,tu_ind);
close(figure(3))
close(figure(7))
close(figure(4))
close(figure(5))
close(figure(9))
close(figure(8))

%% Erk
figure(6)
set(gca,'TickDir','out')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2 1.5])
set(gcf,'renderer','opengl');
feval('print', '-dpdf','-r1000','./project/figures/datafitE.pdf');

% TrkA
figure(2)
set(gca,'TickDir','out')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2 1.5])
set(gcf,'renderer','opengl');
feval('print', '-dpdf','-r1000','./project/figures/datafitT.pdf');

% pErk
figure(1)
for i = 3:-1:1
    s = subplot(1,3,i);
    set(gca,'ylim',[0,0.15])
    
end
set(gca,'TickDir','out')
ah=axes('position',[.13,.11,0.775,0],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.5,'Color','k');
ah=axes('position',[.13,.11,0,0.8],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.5,'Color','k');
figure(1)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4.5 1.5])
set(gcf,'renderer','opengl');
feval('print', '-dpdf','-r1000','./project/figures/datafitP.pdf');
