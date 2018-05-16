%% Definition of model
M.name = 'test_neg_binomial'; 
M.n_subpop = 1; 
M.model = @(T,theta,u)simulate_onestage_1D(T,theta,[]); 
M.theta = @(xi,u)[xi(1);xi(2);xi(3)];
M.dthetadxi= @(xi,u)[1,0,0,0;...
	0,1,0,0;...
	0,0,1,0];

% Subpopulation 1 
% Experiment 1 
s=1; e=1;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [];
M.w_ind{s,e} = [];
M.tau{s,e} = @(t,x,rho,xi,u) [-(x(:,1).*(rho - 1))/rho];
M.dtaudxi{s,e} = @(t,x,dxdxi,rho,drhodxi,xi,u) bsxfun(@times,rho.^(-2), bsxfun(@times,permute(dxdxi(1,:,:),[3,2,1]),(1-rho))-bsxfun(@times,x(:,1),drhodxi));

M.rho{s,e} = @(t,x,xi,u)	[xi(4);xi(4);xi(4);xi(4);xi(4)];
M.drhodxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@plus, [0;0;0;0;0],...
	[0,0,0,1;...
	0,0,0,1;...
	0,0,0,1;...
	0,0,0,1;...
	0,0,0,1]);

M.w{s,e} = @(t,x,xi,u)[1;1;1;1;1];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0;0],...
	[0,0,0,0;...
	0,0,0,0;...
	0,0,0,0;...
	0,0,0,0;...
	0,0,0,0]);
M.distribution{s,e} = 'neg_binomial';

M.u{s,e} = [1];
r=1;
M.scaling{r,e} = @(xi,u)[1];
M.offset{r,e} = @(xi,u)[0];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0];


parameters.name = {'log_{10}(k_{1})',...
'log_{10}(k_{2})',...
'log_{10}(k_3)',...
'\rho'};
parameters.number = length(parameters.name);
parameters.max = [3;3;3;1];
parameters.min = [-3;-3;-3;0];
