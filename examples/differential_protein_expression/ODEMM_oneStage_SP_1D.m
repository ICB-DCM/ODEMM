%% Definition of model
M.name = 'oneStage_SP_1D'; 
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
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [3];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(8))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(8)),[0,0,0,0,0,0,0,2*10^(2*xi(8))*log(10),0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(8))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[xi(9);xi(9);xi(9);xi(9)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0],...
	[0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,1]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [1];
r=1;
M.scaling{r,e} = @(xi,u)[1];
M.offset{r,e} = @(xi,u)[0];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0];
% Subpopulation 1 
% Experiment 2 
s=1; e=2;
M.mean_ind{s,e} = [2];
M.var_ind{s,e} = [5];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(8))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(8)),[0,0,0,0,0,0,0,2*10^(2*xi(8))*log(10),0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(8))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[xi(9);xi(9);xi(9);xi(9)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0],...
	[0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,1]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [1];
r=1;
M.scaling{r,e} = @(xi,u)[1];
M.offset{r,e} = @(xi,u)[0];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 1 
s=2; e=1;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [3];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(8))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(8)),[0,0,0,0,0,0,0,2*10^(2*xi(8))*log(10),0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(8))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[1 - xi(9);1 - xi(9);1 - xi(9);1 - xi(9)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0],...
	[0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,-1]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [0];
r=1;
M.scaling{r,e} = @(xi,u)[1];
M.offset{r,e} = @(xi,u)[0];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 2 
s=2; e=2;
M.mean_ind{s,e} = [2];
M.var_ind{s,e} = [5];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(8))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(8)),[0,0,0,0,0,0,0,2*10^(2*xi(8))*log(10),0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(8))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[1 - xi(9);1 - xi(9);1 - xi(9);1 - xi(9)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0],...
	[0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,-1]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [0];
r=1;
M.scaling{r,e} = @(xi,u)[1];
M.offset{r,e} = @(xi,u)[0];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0];


parameters.name = {'log_{10}(\lambda_0)',...
'log_{10}(\lambda_{A,1})',...
'log_{10}(\lambda_{A,2})',...
'log_{10}(\lambda_{B,1})',...
'log_{10}(\lambda_{B,2})',...
'log_{10}(m\gamma)',...
'log_{10}(\Sigma_{\gamma})',...
'log_{10}(\sigma_{noise})',...
'w_1'};
parameters.number = length(parameters.name);
parameters.max = [3;3;3;3;3;3;1;1;1];
parameters.min = [-3;-3;-3;-3;-3;-3;-3;-3;0];
