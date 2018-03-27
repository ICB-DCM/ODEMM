% Functions to set font sizes and colors for the visualization
TextSizes.DefaultAxesFontSize = 6;
TextSizes.DefaultTextFontSize = 6;
set(0,TextSizes);
fs = TextSizes.DefaultTextFontSize;

%% Ex 1 Sigma Points vs Reaction Rate Equations
color.true = [243,144,25]./255;
color.SP = [11,93,24]/255;
color.RRE = [114,115,255]./255;
color.RRE_CI = [0.90,0.90,1;
    0.75,0.75,1;
    0.60,0.60,1;
    0.45,0.45,1;
    0.10,0.10,1];
color.SP_CI = [ 211 237 207
    150,210,138;...
    70,174,69;...
    0, 150, 0
    11,93,24]./255;
color.data = [112,191,65]/255;
%% Ex 2 full vs marginal distribuion
color.CI_marg = color.RRE_CI;
color.CI_full = color.SP_CI;
color.full = color.SP;
color.marg = color.RRE;

%% NGF-induced Erk signaling
color.PDL = color.SP;
color.ColI = [95,50,124]/255;
color.PDL_data = color.data;
color.ColI_data = [155 84 202]/255;

% differences
color.param(1,:) = [0.0244    0.4350    0.8755];
color.param(2,:) = [0.0265    0.6137    0.8135];
color.param(3,:) = [0.1986    0.7214    0.6310];
color.param(4,:) = [0.6473    0.7456    0.4188];
color.param(5,:) = [0.9856    0.7372    0.2537];
color.param(6,:) = [240,225,102]./255;
color.param(7,:) = [0.2081    0.1663    0.5292];

color.Erk = [240,225,102]./255;
color.TrkA = [0.9856    0.7372    0.2537];
color.k1 =  [0.0244    0.4350    0.8755];
color.k2 = [0.0265    0.6137    0.8135];
color.k4 = [0.1986    0.7214    0.6310];
color.k5 =  [0.6473    0.7456    0.4188];

%%
color.MA_1subpop = [245,211,40]./255;
color.MA_1subpop_extr = [179,106,226]./255;

%% optimization and sampling
color.opt = [179,106,226]./255;
color.samp = [245,211,40]./255;
color.samp_edges = [220,189,35]./255;