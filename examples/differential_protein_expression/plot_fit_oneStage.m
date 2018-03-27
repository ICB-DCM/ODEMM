% This script visualizes the model fit for the two-dimensional data.

clear all
close all
clc

%% Visualization Fit
load_plot_settings

load('./results/results_oneStage_SP_2D','M','D','parameters','options');
xi = parameters.MS.par(:,1);

options_plot.x_scale = 'log';
options_plot.data.col{1} = color.data;
options_plot.model.col{1} =  color.full;
options_plot.simulate_musigma = 1;
options_plot.marginals = true;

fh = plotODEMM(D,M,xi,options_plot);

%% Reproduce figure of paper
options_plot.model.level_linewidth = 0.5;
options_plot.model.levelsets{1,1} = [5,15,25];
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.subplot_lin = true;
options_plot.legendflag = false;
options_plot.plainstyle = true;
options_plot.data.kde = true;
options_plot.marginals = false;
options_plot.titleflag = false;

fh = plotODEMM(D,M,xi,options_plot);
%% modify figure for manuscript
figure(fh{1});
mtvals = log10(400:100:2000);
for ind_s = 1:4
    subplot(2,4,ind_s)
    if ind_s > 1
     set(gca, 'ytick', [log10([400:100:2000])],'xtick',...
         log10([400:100:1700]),'xticklabel','',...
         'yticklabel','')
    else       
     set(gca, 'ytick', log10([400:100:2000]),'yticklabel',...
         {'','5','','','','','10','','','','','','','','','',''},...
         'xtick',log10([400:100:1700]),...
         'xticklabel','');
    end
end
s=subplot(2,4,1);
set(s, 'position', [0.1,0.5,0.2,0.25] );
s=subplot(2,4,2);
set(s, 'position', [0.32,0.5,0.2,0.25] );
s=subplot(2,4,3);
set(s, 'position', [0.54,0.5,0.2,0.25] );
s=subplot(2,4,4);
set(s, 'position', [0.76,0.5,0.2,0.25] );
for ind_s = 5:8
    subplot(2,4,ind_s)
     set(gca, 'ytick', [400:100:2000],'yticklabel',{'','5','','','','','10',...
         '','','','','','','','','',''},...
     'xtick', [400:100:1700], 'xticklabel',{'','5','','','','','10','','','','','','',''},...
     'YMinorTick','off','XMinorTick','off')
end
s=subplot(2,4,5);
set(s, 'position', [0.1,0.23,0.2,0.25] );
s=subplot(2,4,6);
set(s, 'position', [0.32,0.23,0.2,0.25] );
set(gca,'yticklabel',{});
s=subplot(2,4,7);
set(s, 'position', [0.54,0.23,0.2,0.25] );
set(gca,'yticklabel',{});
s=subplot(2,4,8);
set(s, 'position', [0.76,0.23,0.2,0.25] );
set(gca,'yticklabel',{});
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.5 4.5])
%feval('print', '-dpdf','-r1000','./figures/datafit_full.pdf');






