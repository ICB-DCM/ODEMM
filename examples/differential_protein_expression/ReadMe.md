# Differential protein expression (one stage)

This example reproduces the results of the analysis of "Identification of differential protein expression using multivariate data" of Loos et al., Cell Systems (2018).

First, the simulation file needs to be compiled using the MATLAB Toolbox AMICI. This can be done by running
`./simulation/generate_simFile_oneStage.m`.

The script `main_oneStage_SP_1D.m` reproduces the results for model calibration bsaed on the marginal distributions. The script `main_oneStage_SP_2D.m` reproduces the results for model calibration bsaed on the marginal distributions. 

The figures of the paper can be reproduced with the scripts `plot_fit_oneStage.m` and `plot_fit_uncertainty.m`.



