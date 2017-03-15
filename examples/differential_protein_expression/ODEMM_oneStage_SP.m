%% Definition of model
M.name = 'oneStage_SP'; 
M.n_subpop = 2; 
M.model = @(T,theta,u)simulate_onestage_SP(T,theta,[]); 
M.theta = @(xi,u)[(2592480341699211*xi(1))/1125899906842624;(2592480341699211*u(2)*xi(2))/1125899906842624 - (2592480341699211*xi(3)*(u(2) - 1))/1125899906842624;(2592480341699211*u(2)*xi(4))/1125899906842624 - (2592480341699211*xi(5)*(u(2) - 1))/1125899906842624;(2592480341699211*xi(6))/1125899906842624;(2592480341699211*xi(7))/562949953421312];
M.dthetadxi= @(xi,u)[2592480341699211/1125899906842624,0,0,0,0,0,0,0,0;...
	0,(2592480341699211*u(2))/1125899906842624,2592480341699211/1125899906842624 - (2592480341699211*u(2))/1125899906842624,0,0,0,0,0,0;...
	0,0,0,(2592480341699211*u(2))/1125899906842624,2592480341699211/1125899906842624 - (2592480341699211*u(2))/1125899906842624,0,0,0,0;...
	0,0,0,0,0,2592480341699211/1125899906842624,0,0,0;...
	0,0,0,0,0,0,2592480341699211/562949953421312,0,0];

% Subpopulation 1 
% Experiment 1 
s=1; e=1;
M.mean_ind{s,e} = [1,2];
M.var_ind{s,e} = [3,4,5];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [log(x(:,1)) - 0.5*Sigma(:,1,1),...
log(x(:,2)) - 0.5*Sigma(:,2,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_logn_mean(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_logn(t,x,xi,2,[10^(2*xi(8));10^(2*xi(8))],'multiplicative');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_logn(t,x,dxdxi,xi,2,[10^(2*xi(8));10^(2*xi(8))],[0,0,0,0,0,0,0,2*10^(2*xi(8))*log(10),0;...
0,0,0,0,0,0,0,2*10^(2*xi(8))*log(10),0],'multiplicative');

M.w{s,e} = @(t,x,xi,u)[xi(9);xi(9);xi(9);xi(9)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0],...
	[0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,1]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [1];
r=1;
M.scaling{r,e} = @(xi,u)[1;1];
M.offset{r,e} = @(xi,u)[0;0];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 1 
s=2; e=1;
M.mean_ind{s,e} = [1,2];
M.var_ind{s,e} = [3,4,5];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [log(x(:,1)) - 0.5*Sigma(:,1,1),...
log(x(:,2)) - 0.5*Sigma(:,2,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_logn_mean(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_logn(t,x,xi,2,[10^(2*xi(8));10^(2*xi(8))],'multiplicative');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_logn(t,x,dxdxi,xi,2,[10^(2*xi(8));10^(2*xi(8))],[0,0,0,0,0,0,0,2*10^(2*xi(8))*log(10),0;...
0,0,0,0,0,0,0,2*10^(2*xi(8))*log(10),0],'multiplicative');

M.w{s,e} = @(t,x,xi,u)[1 - xi(9);1 - xi(9);1 - xi(9);1 - xi(9)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0],...
	[0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,-1]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [0];
r=1;
M.scaling{r,e} = @(xi,u)[1;1];
M.offset{r,e} = @(xi,u)[0;0];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0];


parameters.name = {'log_{10}(\lambda)',...
'log_{10}(\lambda_{A,s1})',...
'log_{10}(\lambda_{A,s2})',...
'log_{10}(\lambda_{B,s1})',...
'log_{10}(\lambda_{B,s2})',...
'log_{10}(m\gamma)',...
'log_{10}(\sigma_{\gamma})',...
'log_{10}(\sigma_{noise})',...
'w_1'};
parameters.number = length(parameters.name);
parameters.max = [3;3;3;3;3;3;1;1;1];
parameters.min = [-3;-3;-3;-3;-3;-3;-3;-3;0];