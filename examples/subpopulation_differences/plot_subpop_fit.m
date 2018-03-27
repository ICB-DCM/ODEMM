% Visualization script for the fits

close all
clear all
clc

load('./results/results_subpop_TrkA','parameters');
load('./data/data_PDL');

load_plot_settings
ODEMM_NGF_subpop_TrkA

options_plot.x_scale = 'log';
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.data.lw = 1;
options_plot.model.lw = 0.8;

for e = 1:4
    options_plot.data.col{e} = color.PDL_data;
    options_plot.model.col{e} =  color.PDL;
end
for d = 1:6
    options_plot.model.levelsets{4,d} = [0.9,1.7,2.5,3.3];
    options_plot.model.levelsets{3,d} = [0.45 0.75 1.15 1.5];
end
options_plot.model.levelsets{4,1} = [1.3,2.6,4,5];
options_plot.model.levelsets{4,3} = [0.9,1.7,2.5,4];

options_plot.boundaries(1).y_min = 1e2; options_plot.boundaries(1).y_max = 2e4;
options_plot.boundaries(2).y_min = 1e2; options_plot.boundaries(2).y_max = 2e4;
options_plot.boundaries(3).y_min = [1e1,1e2]; % 1 x n_meas vector
options_plot.boundaries(3).y_max = [2e4,2e4]; % 1 x n_meas vector
options_plot.boundaries(4).y_min = [9e1,1e2]; % 1 x n_meas vector
options_plot.boundaries(4).y_max = [5e3,2e4]; % 1 x n_meas vector

options_plot.subplot_lin = true;
options_plot.switch_axes = true;
options_plot.marginals = false;
options_plot.simulate_musigma = true;
options_plot.data.markersize = 10;
options_plot.model.level_linewidth = 0.5;
options_plot.data.kde = 1;

xi = parameters.MS.par(:,1);
[conditions,D] = collectConditions(D,M);
fh = plotODEMM(D,M,xi,options_plot);
