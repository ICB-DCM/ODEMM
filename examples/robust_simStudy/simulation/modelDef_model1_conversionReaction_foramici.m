% NOTE: In the manuscript the meaning of k1 and k2 is interchanged!

%% MODEL
% Definition of symbolic variables:
syms A B
syms k1 k2 k3
syms u;
syms time;
syms Omega 
Omega = 100;
% Define state vector:
System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega ];
System.state.variable = [ A; B];
System.state.compartment = { 'cell';'cell'  };
System.state.type     = {'moment';'moment'};
System.state.xmin     = [      0     ;      0     ];
System.state.xmax     = [      inf     ;      inf     ];

System.state.mu0      = [1000-k2*1000/(k3+k2); k2*1000/(k3+k2)];
System.state.C0      = [(500*k2^2 - 500*k2*k3 + 500*k3^2)/(k3*(k2 + k3));...
    -(500*k2^2 - 500*k2*k3 + 500*k3^2)/(k3^2 + k2*k3);...
    (500*k2^2 - 500*k2*k3 + 500*k3^2)/(k3*(k2 + k3))];

System.state.number = 2;
% Define parameter vector:
System.parameter.variable = [ k1;k2;k3];
System.kappa.variable = [];
System.scaleIndicator = 'macroscopic';

% Define propensities:
% A
System.reaction(1).educt      = A;
System.reaction(1).product    = B;
System.reaction(1).propensity = k1*A;

System.reaction(2).educt      = A;
System.reaction(2).product    = B;
System.reaction(2).propensity = k2*A;

System.reaction(3).educt      = B;
System.reaction(3).product    = A;
System.reaction(3).propensity = k3*B;

System.output.variable = [B];
System.output.function = [B];
System.output.name = {'B'};
System.output.length = 1;
