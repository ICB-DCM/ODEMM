clear all
close all
clc

load('./results/results_ECM_differences')
load('./data/data_PDL_ColI'); 

p = parameters{15};
clear parameters;
parameters = p;
ODEMM_NGF_ECM_T_comb15

%% Load and set options
load_plot_settings

options_plot.x_scale = 'log';
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.data.lw = 1;
options_plot.model.lw = 1;

for e = [1,3,5,7]
   options_plot.data.col{e} =  color.PDL_data;
    options_plot.model.col{e} =  color.PDL;
end
for e = [2,4,6,8]
    options_plot.data.col{e} =  color.ColI_data;
    options_plot.model.col{e} = color.ColI;
end
for e = 1:7
    options_plot.model.levelsets{5,e} = [0.45 0.75 1.15 1.5];
    options_plot.model.levelsets{6,e} = [0.45 0.75 1.15 1.5];
    options_plot.model.levelsets{7,e} = [0.9,1.7,2.5,3.3];
    options_plot.model.levelsets{8,e} = [0.9,1.7,2.5,3.3];
end
for e = 1:4
    options_plot.boundaries(e).y_min = 1e2; 
    options_plot.boundaries(e).y_max = 2e4;
end

options_plot.boundaries(5).y_min = [1e0,5e1];
options_plot.boundaries(5).y_max = [1e4,5e5];
options_plot.boundaries(6).y_min = [1e0,5e1];
options_plot.boundaries(6).y_max = [1e4,5e5];
options_plot.boundaries(7).y_min = [1e2,1e2];
options_plot.boundaries(7).y_max = [5e3,2e5];
options_plot.boundaries(8).y_min = [1e2,1e2];
options_plot.boundaries(8).y_max = [5e3,2e5];

options_plot.switch_axes = true;
options_plot.marginals = false;
options_plot.model.subpopulations = false;
options_plot.simulate_musigma = true;
options_plot.data.kde = true;

%% Plot fit
xi = parameters.MS.par(:,1);
[conditions,D] = collectConditions(D,M);
fh = plotODEMM(D,M,xi,options_plot);
