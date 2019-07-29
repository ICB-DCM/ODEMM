init = [];

%% MODEL
% Definition of symbolic variables:
syms A 
syms lambda0 lambdaA gammaAB
syms u;
syms time;

% Define state vector:
System.time = time;
System.compartments = {'cell'};
System.volumes = [1000];
System.state.variable = [ A];
System.state.compartment = { 'cell';  };
System.state.type     = {'moment'};
System.state.xmin     = [      0        ];
System.state.xmax     = [      inf   ];

if isempty(init)
    System.state.mu0      = [    lambda0/gammaAB ];
    System.state.C0      = [ lambda0/gammaAB ];
else
    System.state.mu0      = [init];
end
System.state.number = 2;
% Define parameter vector:
System.parameter.variable = [ lambda0 ;lambdaA;  gammaAB];
System.kappa.variable = [];
System.scaleIndicator = 'microscopic';


% Define propensities:
System.reaction(1).educt      = [];
System.reaction(1).product    = A;
System.reaction(1).propensity = lambda0;

System.reaction(2).educt      = [];
System.reaction(2).product    = A;
System.reaction(2).propensity = lambdaA;

System.reaction(3).educt      = A;
System.reaction(3).product    = [];
System.reaction(3).propensity = gammaAB*A;

System.output.variable = [A];
System.output.function = [A];
System.output.name = {'A'};
System.output.length = 2;
