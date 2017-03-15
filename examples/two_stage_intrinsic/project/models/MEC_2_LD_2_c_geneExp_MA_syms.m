% function [model] = MEC_2_LD_2_c_geneExp_MA_syms(f0_user)
function [model] = MEC_2_LD_2_c_geneExp_MA_syms(varargin)

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;

% STATES

syms mu_1 mu_2 C_1_1 C_1_2 C_2_2

x = [
mu_1, mu_2, C_1_1, C_1_2, C_2_2 ...
];

% PARAMETERS

syms k1 k2 k3 k4 k5 

% KAPPA (constant parameters)

syms indmu1 indmu2 indC1 indC2 indC3 kmu01 kmu02 kC01 kC02 kC03 

syms t

p = [k1,k2,k3,k4,k5];

k = [indmu1,indmu2,indC1,indC2,indC3,kmu01,kmu02,kC01,kC02,kC03];

if nargin > 0
   f0_user = varargin{1};
   if ~isnumeric(f0_user)
      p_user = setdiff(symvar(f0_user),p);
      % ADDITIONAL PARAMETERS IN INITIAL CONDITIONS
      p = [p,p_user];
   end
	fmu01 = f0_user(1); 
	fmu02 = f0_user(2); 
	fC01 = f0_user(3); 
	fC02 = f0_user(4); 
	fC03 = f0_user(5); 
else
	fmu01 = k1/k3; 
	fmu02 = (k1*k4)/(k3*k5); 
	fC01 = k1/k3; 
	fC02 = (k1*k4)/(k3*(k3 + k5)); 
	fC03 = (k1*k4)/(k3*k5) + (k1*k4^2)/(k3*k5*(k3 + k5)); 
end
% INPUT 

u = sym.empty(0,0);

% SYSTEM EQUATIONS

xdot = sym(zeros(size(x)));

xdot(1) = k1/1000 + k2/1000 - k3*mu_1;
xdot(2) = k4*mu_1 - k5*mu_2;
xdot(3) = k1/1000000 + k2/1000000 - 2*C_1_1*k3 + (k3*mu_1)/1000;
xdot(4) = C_1_1*k4 - C_1_2*k3 - C_1_2*k5;
xdot(5) = 2*C_1_2*k4 - 2*C_2_2*k5 + (k4*mu_1)/1000 + (k5*mu_2)/1000;
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = (indmu1*kmu01)/1000 - (fmu01*(indmu1 - 1))/1000;
x0(2) = (indmu2*kmu02)/1000 - (fmu02*(indmu2 - 1))/1000;
x0(3) = (indC1*kC01)/1000000 - (fC01*(indC1 - 1))/1000000;
x0(4) = (indC2*kC02)/1000000 - (fC02*(indC2 - 1))/1000000;
x0(5) = (indC3*kC03)/1000000 - (fC03*(indC3 - 1))/1000000;

% OBSERVABLES

y = sym(zeros(2,1));

y(1) = mu_2;
y(2) = C_2_2;

% SYSTEM STRUCT

model.sym.nmx = 0;
model.sym.x = x;
model.sym.u = u;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.y = y;
% Additional fields for the prespecified length of kappa
model.sym.nk1 = 0;
end