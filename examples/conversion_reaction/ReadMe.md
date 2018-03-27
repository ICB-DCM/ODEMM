# Conversion Reaction

This example reproduces the results of the analysis of "Unraveling sources of heterogeneity" of Loos et al., Cell Systems (2018).

First, the simulation file needs to be compiled using the MATLAB Toolbox AMICI. This can be done by running
`./simulation/generate_simFiles_cr.m`.

All models for this example are optimized within the function `main_optimization_cr.m`.

The Bayesian model selection is performed in `main_sampling_cr.m`.

The prediction of single-cell trajectories is performed in `main_singlecell_prediction.m`.

To reproduce the figures of the paper, `run plot_cr.m`.
