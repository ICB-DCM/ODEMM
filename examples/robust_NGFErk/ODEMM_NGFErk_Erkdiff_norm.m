%% Definition of model
M.name = 'NGFErk_Erkdiff_norm'; 
M.n_subpop = 2; 
M.model = @(T,theta,u)simulate_SigmaPoints_sPsET_loglog_corr(T,theta,u(1)); 
M.theta = @(xi,u)[(2592480341699211*xi(1))/1125899906842624;(2592480341699211*xi(2))/1125899906842624;(2592480341699211*xi(3))/1125899906842624;(2592480341699211*xi(4))/1125899906842624;(2592480341699211*u(2)*xi(5))/1125899906842624 - (2592480341699211*xi(6)*(u(2) - 1))/1125899906842624;(2592480341699211*u(2)*xi(7))/1125899906842624 - (2592480341699211*xi(8)*(u(2) - 1))/1125899906842624;log(log(10^(2*u(2)*xi(9) - 2*xi(10)*(u(2) - 1)) + 1));(2592480341699211*xi(12))/1125899906842624;log(log(10^(2*xi(11)) + 1))];
M.dthetadxi= @(xi,u)[2592480341699211/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,2592480341699211/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,2592480341699211/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,2592480341699211/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,(2592480341699211*u(2))/1125899906842624,2592480341699211/1125899906842624 - (2592480341699211*u(2))/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,(2592480341699211*u(2))/1125899906842624,2592480341699211/1125899906842624 - (2592480341699211*u(2))/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,(2*10^(2*u(2)*xi(9) - 2*xi(10)*(u(2) - 1))*u(2)*log(10))/(log(10^(2*u(2)*xi(9) - 2*xi(10)*(u(2) - 1)) + 1)*(10^(2*u(2)*xi(9) - 2*xi(10)*(u(2) - 1)) + 1)),-(10^(2*u(2)*xi(9) - 2*xi(10)*(u(2) - 1))*log(10)*(2*u(2) - 2))/(log(10^(2*u(2)*xi(9) - 2*xi(10)*(u(2) - 1)) + 1)*(10^(2*u(2)*xi(9) - 2*xi(10)*(u(2) - 1)) + 1)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,2592480341699211/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,(2*10^(2*xi(11))*log(10))/(log(10^(2*xi(11)) + 1)*(10^(2*xi(11)) + 1)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

% Subpopulation 1 
% Experiment 1 
s=1; e=1;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [x(:,1)];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [permute(dxdxi(1,:,:),[3,2,1])],...
	bsxfun(@plus,[0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(x(:,2)+10.^(bsxfun(@times,2,xi(24)))).^(bsxfun(@rdivide,1,2))];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_norm(t,x,dxdxi,xi,10^(2*xi(24)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(24))*log(10),0,0,0,0,0,0],'additive'),2*((x(:,2)+10.^(bsxfun(@times,2,xi(24)))).^(bsxfun(@rdivide,1,2))));

M.w{s,e} = @(t,x,xi,u)[xi(30);xi(30);xi(30);xi(30);xi(30);xi(30);xi(30)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0;0;0;0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]);
M.distribution{s,e} = 'norm';

M.u{s,e} = [0];
r=1;
M.scaling{r,e} = @(xi,u)[1];
M.offset{r,e} = @(xi,u)[10^xi(13)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,10^xi(13)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 1 
% Experiment 2 
s=1; e=2;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [x(:,1)];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [permute(dxdxi(1,:,:),[3,2,1])],...
	bsxfun(@plus,[0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(x(:,2)+10.^(bsxfun(@times,2,xi(25)))).^(bsxfun(@rdivide,1,2))];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_norm(t,x,dxdxi,xi,10^(2*xi(25)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(25))*log(10),0,0,0,0,0],'additive'),2*((x(:,2)+10.^(bsxfun(@times,2,xi(25)))).^(bsxfun(@rdivide,1,2))));

M.w{s,e} = @(t,x,xi,u)[xi(30)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]);
M.distribution{s,e} = 'norm';

M.u{s,e} = [0];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(14)];
M.offset{r,e} = @(xi,u)[10^xi(15)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(14)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(15)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 1 
% Experiment 3 
s=1; e=3;
M.mean_ind{s,e} = [3,1];
M.var_ind{s,e} = [9,6,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [x(:,1),x(:,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_norm(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_norm(t,x,xi,2,[10^(2*xi(28));10^(2*xi(26))],'additive');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_norm(t,x,dxdxi,xi,2,[10^(2*xi(28));10^(2*xi(26))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(28))*log(10),0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(26))*log(10),0,0,0,0],'additive');

M.w{s,e} = @(t,x,xi,u)[xi(30)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]);
M.distribution{s,e} = 'norm';

M.u{s,e} = [0];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(20);10^xi(16)];
M.offset{r,e} = @(xi,u)[10^xi(21);10^xi(17)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(20)*log(10),0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(16)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(21)*log(10),0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(17)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 1 
% Experiment 4 
s=1; e=4;
M.mean_ind{s,e} = [2,1];
M.var_ind{s,e} = [7,5,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [x(:,1),x(:,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_norm(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_norm(t,x,xi,2,[10^(2*xi(29));10^(2*xi(27))],'additive');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_norm(t,x,dxdxi,xi,2,[10^(2*xi(29));10^(2*xi(27))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(29))*log(10),0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(27))*log(10),0,0,0],'additive');

M.w{s,e} = @(t,x,xi,u)[xi(30)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]);
M.distribution{s,e} = 'norm';

M.u{s,e} = [0];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(22);10^xi(18)];
M.offset{r,e} = @(xi,u)[10^xi(23);10^xi(19)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(22)*log(10),0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(18)*log(10),0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(23)*log(10),0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(19)*log(10),0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 1 
s=2; e=1;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [x(:,1)];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [permute(dxdxi(1,:,:),[3,2,1])],...
	bsxfun(@plus,[0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(x(:,2)+10.^(bsxfun(@times,2,xi(24)))).^(bsxfun(@rdivide,1,2))];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_norm(t,x,dxdxi,xi,10^(2*xi(24)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(24))*log(10),0,0,0,0,0,0],'additive'),2*((x(:,2)+10.^(bsxfun(@times,2,xi(24)))).^(bsxfun(@rdivide,1,2))));

M.w{s,e} = @(t,x,xi,u)[1 - xi(30);1 - xi(30);1 - xi(30);1 - xi(30);1 - xi(30);1 - xi(30);1 - xi(30)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0;0;0;0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]);
M.distribution{s,e} = 'norm';

M.u{s,e} = [1];
r=1;
M.scaling{r,e} = @(xi,u)[1];
M.offset{r,e} = @(xi,u)[10^xi(13)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,10^xi(13)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 2 
s=2; e=2;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [x(:,1)];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [permute(dxdxi(1,:,:),[3,2,1])],...
	bsxfun(@plus,[0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(x(:,2)+10.^(bsxfun(@times,2,xi(25)))).^(bsxfun(@rdivide,1,2))];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_norm(t,x,dxdxi,xi,10^(2*xi(25)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(25))*log(10),0,0,0,0,0],'additive'),2*((x(:,2)+10.^(bsxfun(@times,2,xi(25)))).^(bsxfun(@rdivide,1,2))));

M.w{s,e} = @(t,x,xi,u)[1 - xi(30)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]);
M.distribution{s,e} = 'norm';

M.u{s,e} = [1];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(14)];
M.offset{r,e} = @(xi,u)[10^xi(15)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(14)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(15)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 3 
s=2; e=3;
M.mean_ind{s,e} = [3,1];
M.var_ind{s,e} = [9,6,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [x(:,1),x(:,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_norm(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_norm(t,x,xi,2,[10^(2*xi(28));10^(2*xi(26))],'additive');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_norm(t,x,dxdxi,xi,2,[10^(2*xi(28));10^(2*xi(26))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(28))*log(10),0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(26))*log(10),0,0,0,0],'additive');

M.w{s,e} = @(t,x,xi,u)[1 - xi(30)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]);
M.distribution{s,e} = 'norm';

M.u{s,e} = [1];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(20);10^xi(16)];
M.offset{r,e} = @(xi,u)[10^xi(21);10^xi(17)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(20)*log(10),0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(16)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(21)*log(10),0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(17)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 4 
s=2; e=4;
M.mean_ind{s,e} = [2,1];
M.var_ind{s,e} = [7,5,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [x(:,1),x(:,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_norm(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_norm(t,x,xi,2,[10^(2*xi(29));10^(2*xi(27))],'additive');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_norm(t,x,dxdxi,xi,2,[10^(2*xi(29));10^(2*xi(27))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(29))*log(10),0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(27))*log(10),0,0,0],'additive');

M.w{s,e} = @(t,x,xi,u)[1 - xi(30)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]);
M.distribution{s,e} = 'norm';

M.u{s,e} = [1];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(22);10^xi(18)];
M.offset{r,e} = @(xi,u)[10^xi(23);10^xi(19)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(22)*log(10),0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(18)*log(10),0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(23)*log(10),0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(19)*log(10),0,0,0,0,0,0,0,0,0,0,0];


parameters.name = {