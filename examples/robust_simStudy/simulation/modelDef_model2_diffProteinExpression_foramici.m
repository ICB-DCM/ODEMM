init = [];

%% MODEL
% Definition of symbolic variables:
syms A B 
syms lambda0 lambdaA lambdaB gammaAB
syms u;
syms time;

% Define state vector:
System.time = time;
System.compartments = {'cell'};
System.volumes = [1000];
System.state.variable = [ A; B];
System.state.compartment = { 'cell';'cell'  };
System.state.type     = {'moment';'moment'};
System.state.xmin     = [      0     ;      0     ];
System.state.xmax     = [      inf     ;      inf     ];

if isempty(init)
    System.state.mu0      = [    lambda0/gammaAB, lambda0/gammaAB ];
    System.state.C0      = [ lambda0/gammaAB, 0,lambda0/gammaAB ];
else
    System.state.mu0      = [init];
end
System.state.number = 2;
% Define parameter vector:
System.parameter.variable = [ lambda0 ;lambdaA; lambdaB; gammaAB];
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

System.reaction(4).educt      = [];
System.reaction(4).product    = B;
System.reaction(4).propensity = lambda0;

System.reaction(5).educt      = [];
System.reaction(5).product    = B;
System.reaction(5).propensity = lambdaB;

System.reaction(6).educt      = B;
System.reaction(6).product    = [];
System.reaction(6).propensity = gammaAB*B;

System.output.variable = [A;B];
System.output.function = [A;B];
System.output.name = {'A','B'};
System.output.length = 2;
