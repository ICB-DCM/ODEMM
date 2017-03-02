% Visualization script for the Figure 5C


close all
clear all
clc

load ./project/data/data_Lys_1D2D_sepscaled
load('./project/results/results_subpop_final','parameters');
load_plot_settings
ODEMM_NGF_subpop_T_final; % TrkA

%% plotting options
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
options_plot.model.levelsets{4,1} = [1.3,2.6,4,5];
options_plot.model.levelsets{4,3} = [0.9,1.7,2.5,4];
options_plot.model.levelsets{4,6} = [0.9,1.7,2.5,3.3];
options_plot.model.levelsets{3,1} = [0.45 0.75 1.15 1.5];
options_plot.model.levelsets{3,3} = [0.45 0.75 1.15 1.5];
options_plot.model.levelsets{3,6} = [0.45 0.75 1.15 1.5];
options_plot.boundaries(1).y_min = 1e2; options_plot.boundaries(1).y_max = 2e4;
options_plot.boundaries(2).y_min = 1e2; options_plot.boundaries(2).y_max = 2e4;
options_plot.boundaries(3).y_min = [1e1,1e2];
options_plot.boundaries(3).y_max = [2e4,2e4];
options_plot.boundaries(4).y_min = [9e1,1e2];
options_plot.boundaries(4).y_max = [5e3,2e4];
options_plot.subplot_lin = 1;
options_plot.switch_axes = true;
options_plot.marginals = 0;
options_plot.simulate_musigma = true;
options_plot.data.markersize = 10;
options_plot.model.level_linewidth = 0.5;
%%
tu_ind{1} = [1:7];
tu_ind{2} = [1:7];
tu_ind{3} = [1,3,6];
tu_ind{4} = [1,3,6];
xi = parameters.MS.par(:,1);
[conditions,D] = collectConditions(D,M);
fh = plotODEMix(D,M,xi,[1:4],options_plot,tu_ind);
%fh = plotODEMix(D,[],[],[1:4],options_plot,tu_ind);

%% pErk kinetc
figure(fh{1});
s = subplot(1,7,1);
set(gca,'TickDir','out')
ah=axes('position',[.13,.11,0.775,0],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.5,'Color','k');
ah=axes('position',[.13,.11,0,0.8],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.5,'Color','k');
for s_ind = 1:7
    s = subplot(1,7,s_ind);
    ylim([0,0.11])
end
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4.5 1.5])
set(gcf,'renderer','opengl');
feval('print', '-dpdf','-r1000','./project/figures/datafit_pErk_kin.pdf');
%feval('print', '-dpdf','-r1000','./project/figures/data_pErk_kin.pdf');

%% pErk dose response
figure(fh{2});
s = subplot(1,7,1);
set(gca,'TickDir','out')
ah=axes('position',[.13,.11,0.775,0],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.5,'Color','k');
ah=axes('position',[.13,.11,0,0.8],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.5,'Color','k');
for s_ind = 1:7
    s = subplot(1,7,s_ind);
    ylim([0,0.12])
end
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4.5 1.5])
set(gcf,'renderer','opengl');

feval('print', '-dpdf','-r1000','./project/figures/datafit_pErk_dr.pdf');
%feval('print', '-dpdf','-r1000','./project/figures/data_pErk_dr.pdf');

%% TrkA/pErk
figure(fh{3});
s=subplot(2,3,1)
set(s, 'position', [0.1,0.5,0.25,0.25] );
set(gca, 'xtick', [2,3,4] ,'xticklabel',{''});
set(gca, 'ytick', [2,3,4] ,'yticklabel',{'10^2','10^3','10^4'});
mtvals(1:8) = log10([2:9]*10^2);
mtvals(9:16) = log10([2:9]*10^3);
hA = gca;
hA.YAxis.MinorTickValues = mtvals;
hA.XAxis.MinorTickValues = mtvals;
set(gca,'YMinorTick','on','XMinorTick','on')
view(180,90);
set(gca,'xdir','reverse');

s=subplot(2,3,2)
set(s, 'position', [0.36,0.5,0.25,0.25] );
set(gca, 'xtick', [2,3,4] ,'xticklabel',{''});
set(gca, 'ytick', [2,3,4] ,'yticklabel',{''});
hA = gca;
hA.YAxis.MinorTickValues = mtvals;
hA.XAxis.MinorTickValues = mtvals;
set(gca,'YMinorTick','on','XMinorTick','on')
view(180,90);
set(gca,'xdir','reverse');
s=subplot(2,3,3)
set(s, 'position', [0.62,0.5,0.25,0.25] );
set(gca, 'xtick', [2,3,4],'xticklabel',{''});
set(gca, 'ytick', [2,3,4] ,'yticklabel',{''});
hA = gca;
hA.YAxis.MinorTickValues = mtvals;
hA.XAxis.MinorTickValues = mtvals;
set(gca,'YMinorTick','on','XMinorTick','on')
view(180,90);
set(gca,'xdir','reverse');
s=subplot(2,3,4)
set(s, 'position', [0.1,0.23,0.25,0.25] );

set(gca, 'ytick', [1e2,1e3,1e4],'yticklabel',{'10^2','10^3','10^4'});
set(gca, 'xtick', [1e2,1e3,1e4] ,'xticklabel',{'10^2','','10^4'});
view(180,90);
set(gca,'xdir','reverse');

s=subplot(2,3,5)
set(s, 'position', [0.36,0.23,0.25,0.25] );
set(gca, 'xtick', [1e2,1e3,1e4] ,'xticklabel',{'10^2','','10^4'});
set(gca, 'ytick', [1e2,1e3,1e4] ,'yticklabel',{''});

view(180,90);
set(gca,'xdir','reverse');
s=subplot(2,3,6)
set(s, 'position', [0.62,0.23,0.25,0.25] );
set(gca, 'xtick', [1e2,1e3,1e4] ,'xticklabel',{'10^2','','10^4'},'XMinorTick','off','YMinorTick','off');
set(gca, 'ytick', [1e2,1e3,1e4] ,'yticklabel',{''});

view(180,90);
set(gca,'xdir','reverse');

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4.5 4.5])
set(gcf,'renderer','opengl');
print('-dpdf','-r500','./project/figures/datafit_TrkApErk')
%% Erk/pErk
figure(fh{4});
s=subplot(2,3,1)
set(s, 'position', [0.1,0.5,0.25,0.25] );
set(gca, 'xtick', [2,3] ,'xticklabel',{''});
set(gca, 'ytick', [2,3,4] ,'yticklabel',{'10^2','10^3','10^4'});
mtvals(1:8) = log10([2:9]*10^2);
mtvals(9:16) = log10([2:9]*10^3);
hA = gca;
hA.YAxis.MinorTickValues = mtvals;
hA.XAxis.MinorTickValues = mtvals;
set(gca,'YMinorTick','on','XMinorTick','on')
view(180,90);
set(gca,'xdir','reverse');

s=subplot(2,3,2)
set(s, 'position', [0.36,0.5,0.25,0.25] );
set(gca, 'xtick', [2,3] ,'xticklabel',{''});
set(gca, 'ytick', [2,3,4] ,'yticklabel',{''});
hA = gca;
hA.YAxis.MinorTickValues = mtvals;
hA.XAxis.MinorTickValues = mtvals;
set(gca,'YMinorTick','on','XMinorTick','on')
view(180,90);
set(gca,'xdir','reverse');

s=subplot(2,3,3)
set(s, 'position', [0.62,0.5,0.25,0.25] );
set(gca, 'xtick', [2,3],'xticklabel',{''});
set(gca, 'ytick', [2,3,4] ,'yticklabel',{''});
hA = gca;
hA.YAxis.MinorTickValues = mtvals;
hA.XAxis.MinorTickValues = mtvals;
set(gca,'YMinorTick','on','XMinorTick','on')
view(180,90);
set(gca,'xdir','reverse');

s=subplot(2,3,4)
set(s, 'position', [0.1,0.23,0.25,0.25] );
set(gca, 'xtick', [1e2,1e3] ,'xticklabel',{'10^2','10^3'});
set(gca, 'ytick', [1e2,1e3,1e4] ,'yticklabel',{'10^2','10^3','10^4'});
view(180,90);
set(gca,'xdir','reverse');

s=subplot(2,3,5)
set(s, 'position', [0.36,0.23,0.25,0.25] );
set(gca, 'xtick', [1e2,1e3] ,'xticklabel',{'10^2','10^3'});
set(gca, 'ytick', [1e2,1e3,1e4] ,'yticklabel',{''});
view(180,90);
set(gca,'xdir','reverse');

s=subplot(2,3,6)
set(s, 'position', [0.62,0.23,0.25,0.25] );
set(gca, 'xtick', [1e2,1e3] ,'xticklabel',{'10^2','10^3'});
set(gca, 'ytick', [1e2,1e3,1e4] ,'yticklabel',{''});
view(180,90);
set(gca,'xdir','reverse');

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4.5 4.5])
set(gcf,'renderer','opengl');
print('-dpdf','-r500','./project/figures/datafit_ErkpErk')
close all

