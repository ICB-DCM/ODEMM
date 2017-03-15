% This script visualizes the model fit for the two-dimensional data.

clear all
close all
clc

%% Visualization Fit
load_plot_settings

load('./project/results/results_diffgeneexp_2D','M','D','parameters','options');
clear options_plot
xi = parameters.MS.par(:,1);

options_plot.x_scale = 'log';
options_plot.data.plot = 'empty';
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.model.level_linewidth = 0.5;
options_plot.model.levelsets{1,1} = [5,15,25];
options_plot.data.col{1} = color.data;
options_plot.model.col{1} =  color.full;
options_plot.subplot_lin = 1;
options_plot.legendflag = 0;
options_plot.plainstyle = 1;
options_plot.data.kde = 1;
options_plot.marginals = 1;

options_plot.simulate_musigma = 1;

fh = plotODEMix(D,M,xi,[],options_plot);

%% modify figure
figure(fh{1});
mtvals = log10(400:100:2000);

% mtvals(1:8) = log10([2:9]*10^2);
% mtvals(9:16) = log10([2:9]*10^3);
% mtvals(17:24) = log10([1100:100:1900]);
for ind_s = 1:4
    subplot(2,4,ind_s)
%     set(gca, 'ytick', log10([400,1000,2000]),'yticklabel',{},...
%     'xtick', log10([400,1000]),'xticklabel',{},...
%     'YMinorTick','on','XMinorTick','on')
%     hA = gca;
%     hA.YAxis.MinorTickValues = mtvals;
%     hA.XAxis.MinorTickValues = mtvals;
    if ind_s > 1
     set(gca, 'ytick', [log10([400:100:2000])],'xtick',log10([400:100:1700]),'xticklabel','',...
         'yticklabel','')
    else       
     set(gca, 'ytick', log10([400:100:2000]),'yticklabel',{'','5','','','','','10',...
         '','','','','','','','','',''},'xtick',log10([400:100:1700]),...
         'xticklabel','');
    end
end
s=subplot(2,4,1);
set(s, 'position', [0.1,0.5,0.2,0.25] );
s=subplot(2,4,2)
set(s, 'position', [0.32,0.5,0.2,0.25] );
s=subplot(2,4,3)
set(s, 'position', [0.54,0.5,0.2,0.25] );
s=subplot(2,4,4)
set(s, 'position', [0.76,0.5,0.2,0.25] );
%%
for ind_s = 5:8
    subplot(2,4,ind_s)
     set(gca, 'ytick', [400:100:2000],'yticklabel',{'','5','','','','','10',...
         '','','','','','','','','',''},...
     'xtick', [400:100:1700], 'xticklabel',{'','5','','','','','10','','','','','','',''},...
     'YMinorTick','off','XMinorTick','off')
%  try
%      set(gca,'XMinorTickLength',0,'YMinorTickLength',0)
%  catch
%  end
end
s=subplot(2,4,5)
set(s, 'position', [0.1,0.23,0.2,0.25] );
%set(gca,'yticklabel',{'400','','','','','','','','','','','','','','1800','',''});
s=subplot(2,4,6)
set(s, 'position', [0.32,0.23,0.2,0.25] );
set(gca,'yticklabel',{});
s=subplot(2,4,7)
set(s, 'position', [0.54,0.23,0.2,0.25] );
set(gca,'yticklabel',{});
s=subplot(2,4,8)
set(s, 'position', [0.76,0.23,0.2,0.25] );
set(gca,'yticklabel',{});
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.5 4.5])
feval('print', '-dpdf','-r1000','./project/figures/datafit_full.pdf');






