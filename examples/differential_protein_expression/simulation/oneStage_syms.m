function [model] = oneStage_syms()
% This function defines the model of the conversion process with
% the logarithm of the output, which are used for the sigma point
% approximation
%
% Return values:
% model: model struct used with amiwrap

model.param = 'log'; 

% STATES
syms A B

x = [
A, B ...
];

% PARAMETERS
syms lambda_0 lambda_A lambda_B gamma_A gamma_B 

syms t

p = [lambda_0,lambda_A,lambda_B,gamma_A,gamma_B];

k = [];
u = sym.empty(0,0);

% SYSTEM EQUATIONS

xdot = sym(zeros(size(x)));

xdot(1) = lambda_0 + lambda_A - A*gamma_A;
xdot(2) = lambda_0 + lambda_B - B*gamma_B;
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = lambda_0/gamma_A;
x0(2) = lambda_0/gamma_B;

% OBSERVABLES
y = sym(zeros(2,1));

y(1) = log(A);
y(2) = log(B);

% SYSTEM STRUCT
model.sym.x = x;
model.sym.u = u;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.y = y;
end