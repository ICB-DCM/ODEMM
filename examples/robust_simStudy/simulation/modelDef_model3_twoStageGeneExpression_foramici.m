
% if nargin >= 1
%     init = varargin{1};
% else
    init = [];
%end

%% MODEL
% Definition of symbolic variables:
syms mA A
syms k1 k2 k3 k4 k5
syms u;
syms time;

% Define state vector:
System.time = time;
System.compartments = {'cell'};
System.volumes = [1000];
System.state.variable = [ mA; A];
System.state.compartment = { 'cell';'cell'  };
System.state.type     = {'moment';'moment'};
System.state.xmin     = [      0     ;      0     ];
System.state.xmax     = [      inf     ;      inf     ];

if isempty(init)
    System.state.mu0      = [      k1/k3     ;  (k1*k4)/(k3*k5)];
    System.state.C0      = [  k1/k3; k1*k4/(k3*(k3+k5)) ;k4*k1/(k3*k5*(k3+k5)) + (k4*k1)/(k3*k5)];
else
    System.state.mu0      = [init];
end
System.state.number = 2;
% Define parameter vector:
System.parameter.variable = [ k1;k2;k3;k4;k5];
System.kappa.variable = [];
System.scaleIndicator = 'microscopic';


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
