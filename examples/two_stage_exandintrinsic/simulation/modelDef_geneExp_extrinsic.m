% Model definition of the reaction network for the two stage gene
% expression with additional extrinsic noise.

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
% Definition of symbolic variables:
syms mA A
syms k1 k2 k3 k4 k5;
syms mk1 mk2 mk3 mk4 mk5;
syms vk1 vk2 vk3 vk4 vk5;
syms u;
syms time;

% Define state vector:
System.time = time;
System.compartments = {'cell'};
System.volumes = [1000];
System.state.variable = [ mA; A; k1; k2; k3; k4; k5];
System.state.compartment = { 'cell';'cell';'cell';'cell';'cell';'cell';'cell'  };
%System.state.name = {'mA';'A';'mB';'B'}
System.state.type     = {'moment';'moment';'moment';'moment';'moment';'moment';'moment'};
System.state.xmin     = [      0     ;      0   ;0;0;0;0;0  ];
System.state.xmax     = [      inf     ;      inf   ; inf;   inf; inf; inf; inf; ];
syms covmAk1 covmAk2 covmAk3 covmAk4 covmAk5
syms covAk1 covAk2 covAk3 covmA4 covAk5
syms mA0 varA varmA

covmAk1 = vk1/mk3; 
covmAk2 = 0;
covmAk4 = 0;
covmAk5 = 0;


mA0 = mk1*mk3/(mk3^2-vk3);%-covmAk3*mk3/(vk3);
covmAk3= -vk3*mA0/mk3;%mk1*vk3/(vk3-mk3^2);

A0 = (mA0*mk4/mk5)/(1-vk5/mk5^2);%mA0*mk4/(vk5/mk5+mk5);

covAk1 = mk4*covmAk1/mk5;
covAk2 = 0;
covAk3 = mk4*covmAk3/mk5;
covAk4 = vk4*mA0/mk5;
covAk5 = -vk5*A0/mk5;


varmA = (2*covmAk1+covmAk3/1000 + mk1/1000 - 2*covmAk3*mA0+ mA0*mk3/1000)/(2*mk3);

covmAA = (covAk1+mk4*varmA + covmAk4*mA0- covAk3*mA0)/(mk3+mk5);

varA = (covmAk4/1000+covAk5/1000+2*covmAA*mk4 + 2*covAk4*mA0 - 2*covAk5*A0+mA0*mk4/1000+A0*mk5/1000)/(2*mk5);


System.state.mu0      = [mA0*1000;...
    A0*1000;...
    mk1;...
    mk2;...
    mk3;...
    mk4;...
    mk5;...
    ];
System.state.C0      = [ varmA*1e6; %11
    covmAA*1e6; %12
    covmAk1*1000;%13
    covmAk2*1000;%14
    covmAk3*1000;%15
    covmAk4*1000;%16
    covmAk5*1000;%17
    varA*1e6; %22
    covAk1*1000; %23
    covAk2*1000; %24
    covAk3*1000; %25
    covAk4*1000; %26
    covAk5*1000; %27
    vk1; %33
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    vk3;
    0;
    0;
    vk4;
    0;
    vk5];


System.state.number = 7;
% Define parameter vector:
System.parameter.variable = [ mk1;mk2;mk3;mk4;mk5; vk1;vk2;vk3;vk4;vk5];
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
