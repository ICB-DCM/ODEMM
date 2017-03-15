# ODEMMToolbox - ODE constrained mixture modeling

ODEMMToolbox is a toolbox for modeling heterogeneous cell populations in MathWorks MATLAB using mixture models and ordinary differential equation (ODE) models. 

ODEMMToolbox offers 
* models that are able to include mechanistic descriptions of the means of individual subpopulations, and 
* hierarchical population models that incorporate means and covariances of the individual subpopulations. 

These function are demonstrated in several examples included in the [`examples/`](examples/) directory.

## Installation

If the zip archive was downloaded, it needs to be unzipped and the main folder has to be added to the MATLAB search path (non-recursively). 

If the repository was cloned, the main folder needs to be added to the MATLAB search path (non-recursively).

*Note:* Detailed instructions on how to modify your MATLAB search path are provided here: https://de.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html

## Licensing

See [ ```LICENSE```](LICENSE) file in the ODEMM source directory.

## Documentation

ODEMMToolbox usage is demonstrated in various [examples](examples/). Further documentation is available in [```doc/ODEMM-doc.pdf```](doc/ODEMM-doc.pdf).


## Examples 

The following examples are included: 
* Conversion reaction (`examples/conversion_reaction/`)
* Differential protein expression (`examples/differential_protein_expression`)
* Subpopulation differences in NGF-induced Erk1/2 signaling (`examples/subpopulation_differences`)
* Differences mediated by extracellular scaffolds in NGF-induced Erk1/2 signaling (`examples/ECM_differences`)
* Intrinsic noise for a two stage model of gene expression (`examples/two_stage_intrinsic`)


These models require the freely available toolboxes 
* AMICI for simulation (http://icb-dcm.github.io/AMICI/) 
* PESTO for the parameter estimation (http://icb-dcm.github.io/PESTO/)

* SPToolbox for the sigma-point approximation (http://icb-dcm.github.io/SPToolbox/)
