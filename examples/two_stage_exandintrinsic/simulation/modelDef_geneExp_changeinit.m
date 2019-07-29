% Model definition of the reaction network
% This routine should return a struct, called System, with the following fields:
%       System.time 
%       System.compartments
%       System.volumes
%       System.state.variable
%       System.state.compartment
%       System.state.number   = length(System.state.variable);
%       System.state.type
%       System.state.name
%       System.state.xmin
%       System.state.xmax
%       System.state.mu0 
%       System.state.C0  
%       System.state.constraint 
%       System.parameter.variable 
%       System.parameter.MacroscopicValue
%       System.parameter.MicroscopicValue
%       System.parameter.name 
%       System.reaction.educt  
%       System.reaction.product
%       System.reaction.propensity  : Microscopic definition of the
%       propensity
%       System.reaction.MacroscopicPropensity: OPTIONAL. Macroscopic
%       definition of the propensity.

%% MODEL
function System = modelDef_geneExp_changeinit(init)

% Definition of symbolic variables:
syms mA A 
syms k1 k2 k3 k4 k5;
syms u;
syms time;

% Define state vector:
System.time = time;
System.compartments = {'cell'};
System.volumes = [1000];
System.state.variable = [ mA; A];
System.state.compartment = { 'cell';'cell'  };
%System.state.name = {'mA';'A';'mB';'B'}
System.state.type     = {'moment';'moment'};
System.state.xmin     = [      0     ;      0     ];
System.state.xmax     = [      inf     ;      inf     ];
System.state.mu0      = [init];
System.state.number = 2;
% Define parameter vector:
System.parameter.variable = [ k1;k2;k3;k4;k5];
System.kappa.variable = [];
System.scaleIndicator = 'macroscopic';


% Define propensities:
% A
System.reaction(1).educt      = [];
System.reaction(1).product    = mA;
System.reaction(1).propensity = k1;

System.reaction(2).educt      = [];
System.reaction(2).product    = mA;
System.reaction(2).propensity = k2;

System.reaction(3).educt      = mA;
System.reaction(3).product    = [];
System.reaction(3).propensity = k3*mA;

System.reaction(4).educt      = mA;
System.reaction(4).product    = [mA;A];
System.reaction(4).propensity = k4*mA;

System.reaction(5).educt      = A;
System.reaction(5).product    = [];
System.reaction(5).propensity = k5*A;

System.output.variable = [A];
System.output.function = [A];
System.output.name = {'A'};
System.output.length = 1;
