function model = ODEmodel_sPsET_loglog_syms()
% This function provides the model of NGF-induced Erk1/2 signaling. It is
% compiled using AMICI.

model.param = 'log';

%% States

syms x1 x2

x = [x1 x2]; %k3*TrkA0:NGF, spErk0

syms p1 p2 p3 p4 p5 p6 %k1 k2 k4 k5 k3TrkA0 sErk0
p = [p1 p2 p3 p4 p5 p6];

syms u1

k = u1;

xdot = sym(zeros(size(x)));

xdot(1) = p1*u1*(p5-x1)-p2*x1;
xdot(2) = (p3 + x1)*(p6-x2)-p4*x2;

x0 = sym(zeros(size(x)));

x0(1) = 0;
x0(2) = (p3*p6)/(p3+p4);

y = sym(zeros(1,1));
y = log([x2,p6,p5]);

model.sym.x = x;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.y = y;

end