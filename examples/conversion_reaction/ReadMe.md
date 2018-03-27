# Conversion Reaction

This example reproduces the results of the analysis of "Unraveling sources of heterogeneity" of Loos et al., Cell Systems (2018). In this example, models incorporating only the **mean** (model_RRE.m), and models incorporating **mean and variance** (e.g., model_SP_k3.m) are implemented. The latter allow for cell-to-cell variability of certain parameters of the model and incorporate the sigma-point approximation to obtain the statistical properties of the individual subpopulations.

![](../images/cr.jpg "Model for the conversion process.")
@image latex "cr.jpg" "Model for the conversion process." width=0.7\textwidth

First, the simulation file needs to be compiled using the MATLAB Toolbox AMICI. This can be done by running
`./simulation/generate_simFiles_cr.m`.

All models for this example are optimized within the function `main_optimization_cr.m`.

The Bayesian model selection is performed in `main_sampling_cr.m`.

The prediction of single-cell trajectories is performed in `main_singlecell_prediction.m`.

To reproduce the figures of the paper, `run_plot_cr.m`.


