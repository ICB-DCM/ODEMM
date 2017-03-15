% This script generates the figures for the supplement

clear all
close all
clc
load_plot_settings
%% Figure SI Appendix
options_plot.x_scale = 'lin';
options_plot.data.plot = 'filled';
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.boundaries(1).x_min = [0.4]; 
options_plot.boundaries(1).x_max = [0.8];
options_plot.subplot_lin = 1;
options_plot.plainstyle =  0;
options_plot.legendflag =  1;
options_plot.model.lw =  1;
options_plot.data.lw =  0.7;
options_plot.data.fill_col{1} = color.data;
options_plot.model.col{1} = color.SP;
options_plot.boundaries(1).y_min = 0.35;
options_plot.boundaries(1).y_max = 0.83;

models = {'RRE_onlyone','RRE_subpop','RRE_timedep','SP_k1', 'SP_k1k2', 'SP_all', 'SP_k2', 'SP_k2k3', 'SP_k1k3', 'SP_k3'};
for m = 1:10
    load(['./project/results/results_' models{m}],'parameters','M','D')
    if m > 3
        options_plot.simulate_musigma = 1; % for Sigma Points
    else
        options_plot.simulate_musigma = 0; % for RRE
    end   
    xi = parameters.MS.par(:,1);    
    plotODEMix(D,M,xi,[],options_plot);
    for i = 1:5
        s = subplot(1,5,i);
        ylim([0,0.15])
    end
    subplot(1,5,1);
    set(gca,'xtick',[0.4,0.8],...
        'xticklabel',{'0.4','0.8'},...
        'tickdir','out');
    ah=axes('position',[.13,.11,0.7813,0],'visible','off');
    line([0,1],[0,1],'parent',ah,'linewidth',0.6,'Color','k');
    ah=axes('position',[.133,.11,0,0.85],'visible','off');
    line([0,1],[0,1],'parent',ah,'linewidth',0.6,'Color','k');
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4 2])
    set(gcf,'renderer','opengl');
    feval('print', '-depsc','-r1000',['./project/figures/datafit_cr_' num2str(m) '.eps']);
end
