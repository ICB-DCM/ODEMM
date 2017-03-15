% Model file for combination 15, i.e. differences in TrkA activity, Erk1/2
% levels and Erk1/2 inactivation

%% Definition of model
M.name = 'NGF_ECM_T_comb15'; 
M.n_subpop = 2; 
M.model = @(T,theta,u)simulate_SigmaPoints_sPsET_loglog_corr(T,theta,u(1)); 
M.theta = @(xi,u)[(2592480341699211*xi(1))/1125899906842624;(2592480341699211*xi(2))/1125899906842624;(2592480341699211*xi(3))/1125899906842624;(2592480341699211*xi(4))/1125899906842624 + (2592480341699211*u(3)*xi(30))/1125899906842624;(2592480341699211*u(2)*xi(5))/1125899906842624 + (2592480341699211*u(3)*xi(31))/1125899906842624 - (2592480341699211*xi(6)*(u(2) - 1))/1125899906842624;(2592480341699211*xi(7))/1125899906842624 + (2592480341699211*u(3)*xi(32))/1125899906842624;log(log(10^(2*u(2)*xi(8) - 2*xi(9)*(u(2) - 1)) + 1));(2592480341699211*xi(11))/1125899906842624;log(log(10^(2*xi(10)) + 1))];
M.dthetadxi= @(xi,u)[2592480341699211/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,2592480341699211/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,2592480341699211/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,2592480341699211/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(2592480341699211*u(3))/1125899906842624,0,0;...
	0,0,0,0,(2592480341699211*u(2))/1125899906842624,2592480341699211/1125899906842624 - (2592480341699211*u(2))/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(2592480341699211*u(3))/1125899906842624,0;...
	0,0,0,0,0,0,2592480341699211/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(2592480341699211*u(3))/1125899906842624;...
	0,0,0,0,0,0,0,(2*10^(2*u(2)*xi(8) - 2*xi(9)*(u(2) - 1))*u(2)*log(10))/(log(10^(2*u(2)*xi(8) - 2*xi(9)*(u(2) - 1)) + 1)*(10^(2*u(2)*xi(8) - 2*xi(9)*(u(2) - 1)) + 1)),-(10^(2*u(2)*xi(8) - 2*xi(9)*(u(2) - 1))*log(10)*(2*u(2) - 2))/(log(10^(2*u(2)*xi(8) - 2*xi(9)*(u(2) - 1)) + 1)*(10^(2*u(2)*xi(8) - 2*xi(9)*(u(2) - 1)) + 1)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,2592480341699211/1125899906842624,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,(2*10^(2*xi(10))*log(10))/(log(10^(2*xi(10)) + 1)*(10^(2*xi(10)) + 1)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

% Subpopulation 1 
% Experiment 1 
s=1; e=1;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(23))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(23)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(23))*log(10),0,0,0,0,0,0,0,0,0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(23))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[10^xi(29);10^xi(29);10^xi(29);10^xi(29);10^xi(29);10^xi(29);10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0;0;0;0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [0;0];
r=1;
M.scaling{r,e} = @(xi,u)[1];
M.offset{r,e} = @(xi,u)[10^xi(12)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,10^xi(12)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 1 
% Experiment 2 
s=1; e=2;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(23))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(23)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(23))*log(10),0,0,0,0,0,0,0,0,0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(23))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[10^xi(29);10^xi(29);10^xi(29);10^xi(29);10^xi(29);10^xi(29);10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0;0;0;0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [0;1];
r=1;
M.scaling{r,e} = @(xi,u)[1];
M.offset{r,e} = @(xi,u)[10^xi(12)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,10^xi(12)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 1 
% Experiment 3 
s=1; e=3;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(24))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(24)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(24))*log(10),0,0,0,0,0,0,0,0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(24))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [0;0];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(17)];
M.offset{r,e} = @(xi,u)[10^xi(18)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(17)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(18)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 1 
% Experiment 4 
s=1; e=4;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(24))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(24)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(24))*log(10),0,0,0,0,0,0,0,0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(24))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [0;1];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(17)];
M.offset{r,e} = @(xi,u)[10^xi(18)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(17)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(18)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 1 
% Experiment 5 
s=1; e=5;
M.mean_ind{s,e} = [3,1];
M.var_ind{s,e} = [9,6,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [log(x(:,1)) - 0.5*Sigma(:,1,1),...
log(x(:,2)) - 0.5*Sigma(:,2,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_logn_mean(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_logn(t,x,xi,2,[10^(2*xi(25));10^(2*xi(26))],'multiplicative');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_logn(t,x,dxdxi,xi,2,[10^(2*xi(25));10^(2*xi(26))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(25))*log(10),0,0,0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(26))*log(10),0,0,0,0,0,0],'multiplicative');

M.w{s,e} = @(t,x,xi,u)[10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [0;0];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(13);10^xi(19)];
M.offset{r,e} = @(xi,u)[10^xi(14);10^xi(20)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,10^xi(13)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(19)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(14)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(20)*log(10),0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 1 
% Experiment 6 
s=1; e=6;
M.mean_ind{s,e} = [3,1];
M.var_ind{s,e} = [9,6,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [log(x(:,1)) - 0.5*Sigma(:,1,1),...
log(x(:,2)) - 0.5*Sigma(:,2,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_logn_mean(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_logn(t,x,xi,2,[10^(2*xi(25));10^(2*xi(26))],'multiplicative');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_logn(t,x,dxdxi,xi,2,[10^(2*xi(25));10^(2*xi(26))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(25))*log(10),0,0,0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(26))*log(10),0,0,0,0,0,0],'multiplicative');

M.w{s,e} = @(t,x,xi,u)[10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [0;1];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(13);10^xi(19)];
M.offset{r,e} = @(xi,u)[10^xi(14);10^xi(20)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,10^xi(13)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(19)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(14)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(20)*log(10),0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 1 
% Experiment 7 
s=1; e=7;
M.mean_ind{s,e} = [2,1];
M.var_ind{s,e} = [7,5,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [log(x(:,1)) - 0.5*Sigma(:,1,1),...
log(x(:,2)) - 0.5*Sigma(:,2,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_logn_mean(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_logn(t,x,xi,2,[10^(2*xi(27));10^(2*xi(28))],'multiplicative');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_logn(t,x,dxdxi,xi,2,[10^(2*xi(27));10^(2*xi(28))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(27))*log(10),0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(28))*log(10),0,0,0,0],'multiplicative');

M.w{s,e} = @(t,x,xi,u)[10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [0;0];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(15);10^xi(21)];
M.offset{r,e} = @(xi,u)[10^xi(16);10^xi(22)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(15)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(21)*log(10),0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(16)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(22)*log(10),0,0,0,0,0,0,0,0,0,0];
% Subpopulation 1 
% Experiment 8 
s=1; e=8;
M.mean_ind{s,e} = [2,1];
M.var_ind{s,e} = [7,5,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [log(x(:,1)) - 0.5*Sigma(:,1,1),...
log(x(:,2)) - 0.5*Sigma(:,2,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_logn_mean(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_logn(t,x,xi,2,[10^(2*xi(27));10^(2*xi(28))],'multiplicative');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_logn(t,x,dxdxi,xi,2,[10^(2*xi(27));10^(2*xi(28))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(27))*log(10),0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(28))*log(10),0,0,0,0],'multiplicative');

M.w{s,e} = @(t,x,xi,u)[10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [0;1];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(15);10^xi(21)];
M.offset{r,e} = @(xi,u)[10^xi(16);10^xi(22)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(15)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(21)*log(10),0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(16)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(22)*log(10),0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 1 
s=2; e=1;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(23))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(23)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(23))*log(10),0,0,0,0,0,0,0,0,0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(23))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[1 - 10^xi(29);1 - 10^xi(29);1 - 10^xi(29);1 - 10^xi(29);1 - 10^xi(29);1 - 10^xi(29);1 - 10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0;0;0;0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [1;0];
r=1;
M.scaling{r,e} = @(xi,u)[1];
M.offset{r,e} = @(xi,u)[10^xi(12)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,10^xi(12)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 2 
s=2; e=2;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(23))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(23)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(23))*log(10),0,0,0,0,0,0,0,0,0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(23))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[1 - 10^xi(29);1 - 10^xi(29);1 - 10^xi(29);1 - 10^xi(29);1 - 10^xi(29);1 - 10^xi(29);1 - 10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0;0;0;0;0;0;0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [1;1];
r=1;
M.scaling{r,e} = @(xi,u)[1];
M.offset{r,e} = @(xi,u)[10^xi(12)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,10^xi(12)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 3 
s=2; e=3;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(24))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(24)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(24))*log(10),0,0,0,0,0,0,0,0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(24))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[1 - 10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [1;0];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(17)];
M.offset{r,e} = @(xi,u)[10^xi(18)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(17)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(18)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 4 
s=2; e=4;
M.mean_ind{s,e} = [1];
M.var_ind{s,e} = [4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,sigma,xi,u) [log(x(:,1)) - sigma.^2/2];
M.dmudxi{s,e} = @(t,x,dxdxi,sigma,dsigmadxi,xi,u) bsxfun(@plus, [bsxfun(@rdivide, permute(dxdxi(1,:,:),[3,2,1]), x(:,1))],...
	bsxfun(@plus,[bsxfun(@times, -dsigmadxi, sigma)],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]));

M.sigma{s,e} = @(t,x,xi,u)	[(log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(24))).^(1./2)];
M.dsigmadxi{s,e} = @(t,x,dxdxi,xi,u)	bsxfun(@rdivide,func_dsigma2dxi_logn(t,x,dxdxi,xi,10^(2*xi(24)),[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(24))*log(10),0,0,0,0,0,0,0,0],'multiplicative'),2*((log(x(:,2)./x(:,1).^2 + 1) + 10.^(2*xi(24))).^(1./2)));

M.w{s,e} = @(t,x,xi,u)[1 - 10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [1;1];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(17)];
M.offset{r,e} = @(xi,u)[10^xi(18)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(17)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(18)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 5 
s=2; e=5;
M.mean_ind{s,e} = [3,1];
M.var_ind{s,e} = [9,6,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [log(x(:,1)) - 0.5*Sigma(:,1,1),...
log(x(:,2)) - 0.5*Sigma(:,2,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_logn_mean(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_logn(t,x,xi,2,[10^(2*xi(25));10^(2*xi(26))],'multiplicative');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_logn(t,x,dxdxi,xi,2,[10^(2*xi(25));10^(2*xi(26))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(25))*log(10),0,0,0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(26))*log(10),0,0,0,0,0,0],'multiplicative');

M.w{s,e} = @(t,x,xi,u)[1 - 10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [1;0];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(13);10^xi(19)];
M.offset{r,e} = @(xi,u)[10^xi(14);10^xi(20)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,10^xi(13)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(19)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(14)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(20)*log(10),0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 6 
s=2; e=6;
M.mean_ind{s,e} = [3,1];
M.var_ind{s,e} = [9,6,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [log(x(:,1)) - 0.5*Sigma(:,1,1),...
log(x(:,2)) - 0.5*Sigma(:,2,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_logn_mean(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_logn(t,x,xi,2,[10^(2*xi(25));10^(2*xi(26))],'multiplicative');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_logn(t,x,dxdxi,xi,2,[10^(2*xi(25));10^(2*xi(26))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(25))*log(10),0,0,0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(26))*log(10),0,0,0,0,0,0],'multiplicative');

M.w{s,e} = @(t,x,xi,u)[1 - 10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [1;1];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(13);10^xi(19)];
M.offset{r,e} = @(xi,u)[10^xi(14);10^xi(20)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,10^xi(13)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(19)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(14)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(20)*log(10),0,0,0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 7 
s=2; e=7;
M.mean_ind{s,e} = [2,1];
M.var_ind{s,e} = [7,5,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [log(x(:,1)) - 0.5*Sigma(:,1,1),...
log(x(:,2)) - 0.5*Sigma(:,2,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_logn_mean(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_logn(t,x,xi,2,[10^(2*xi(27));10^(2*xi(28))],'multiplicative');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_logn(t,x,dxdxi,xi,2,[10^(2*xi(27));10^(2*xi(28))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(27))*log(10),0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(28))*log(10),0,0,0,0],'multiplicative');

M.w{s,e} = @(t,x,xi,u)[1 - 10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [1;0];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(15);10^xi(21)];
M.offset{r,e} = @(xi,u)[10^xi(16);10^xi(22)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(15)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(21)*log(10),0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(16)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(22)*log(10),0,0,0,0,0,0,0,0,0,0];
% Subpopulation 2 
% Experiment 8 
s=2; e=8;
M.mean_ind{s,e} = [2,1];
M.var_ind{s,e} = [7,5,4];
M.w_ind{s,e} = [];
M.mu{s,e} = @(t,x,Sigma,xi,u) [log(x(:,1)) - 0.5*Sigma(:,1,1),...
log(x(:,2)) - 0.5*Sigma(:,2,2)];

M.dmudxi{s,e} = @(t,x,dxdxi,Sigma,dSigmadxi,xi,u) func_dmudxi_logn_mean(t,x,dxdxi,Sigma,dSigmadxi,xi,u,2);

M.Sigma{s,e} = @(t,x,xi,u) func_Sigma_logn(t,x,xi,2,[10^(2*xi(27));10^(2*xi(28))],'multiplicative');
M.dSigmadxi{s,e} = @(t,x,dxdxi,xi,u) func_dSigmadxi_logn(t,x,dxdxi,xi,2,[10^(2*xi(27));10^(2*xi(28))],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(27))*log(10),0,0,0,0,0;...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2*10^(2*xi(28))*log(10),0,0,0,0],'multiplicative');

M.w{s,e} = @(t,x,xi,u)[1 - 10^xi(29)];
M.dwdxi{s,e} = @(t,x,dxdxi,xi,u)bsxfun(@plus, [0],...
	[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10^xi(29)*log(10),0,0,0]);
M.distribution{s,e} = 'logn_mean';

M.u{s,e} = [1;1];
r=1;
M.scaling{r,e} = @(xi,u)[10^xi(15);10^xi(21)];
M.offset{r,e} = @(xi,u)[10^xi(16);10^xi(22)];
M.dscalingdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(15)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(21)*log(10),0,0,0,0,0,0,0,0,0,0,0];
M.doffsetdxi{r,e} = @(xi,u)[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(16)*log(10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10^xi(22)*log(10),0,0,0,0,0,0,0,0,0,0];


parameters.name = {'log_{10}(k_1)',...
'log_{10}(k_2)',...
'log_{10}(k_4)',...
'log_{10}(k_{5})',...
'log_{10}(k_3TrkA0_{s1})',...
'log_{10}(k_3TrkA0_{s2})',...
'log_{10}(cP_1*mErk0_{s1})',...
'log_{10}(cvTrkA0_s1)',...
'log_{10}(cvTrkA0_s2)',...
'log_{10}(cvErk_1)',...
'log_{10}(corr)',...
'log_{10}(cP_1)',...
'log_{10}(cT_5/k3)',...
'log_{10}(oT_5)',...
'log_{10}(cE_7/cP_5})',...
'log_{10}(bE_7)',...
'log_{10}(cP_3)',...
'log_{10}(oP_3)',...
'log_{10}(cP_5)',...
'log_{10}(oP_5)',...
'log_{10}(cP_7)',...
'log_{10}(oP_7)',...
'log_{10}(\epsilon_P_{12})',...
'log_{10}(\epsilon_P_{34})',...
'log_{10}(\epsilon_T_{56})',...
'log_{10}(\epsilon_P_{56})',...
'log_{10}(\epsilon_E_{78})',...
'log_{10}(\epsilon_P_{78})',...
'log_{10}(w)',...
'log_{10}(\kappa_{k_5})',...
'log_{10}(\kappa_{k_3TrkA_0})',...
'log_{10}(\kappa_{c_{P1}Erk})'};
parameters.number = length(parameters.name);
parameters.max = [6;6;6;6;6;6;6;4;4;4;4;6;6;6;6;6;6;6;6;6;6;6;1;1;1;1;1;1;0;3;3;3];
parameters.min = [-6;-6;-6;-6;-6;-6;-6;-4;-4;-4;-4;-6;-6;-6;-6;-6;-6;-6;-6;-6;-6;-6;-3;-3;-3;-3;-3;-3;-4;-3;-3;-3];