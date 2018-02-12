close all
clear all
clc

xi_true = [1.7,2.7,2,2,2.7,  -1,-1,-1,0.5];
% lambda_0, lambda_A,s1,lambda_A,s2, lambda_B,s1, lambda_B,s2, m_gamma, sigma_gamma, sigma_noise, w                   

w = 0.5;
sigma_gamma = 10.^(xi_true(7));


t = [0,0.5,1,2]; % hours

n_data = 1000;
theta_s1 = log(10)*xi_true([1,2,4,6,6]); % same mean degradation in both subpopulations, log parameterized for SP
theta_s2 = log(10)*xi_true([1,3,5,6,6]);
m_gamma = theta_s1(4);
sigma_noise = 10^(xi_true(end-1));
D(1).y = nan(1,numel(t),n_data,2);
D(1).t = t;
D(1).n_dim = 2;
D(1).measurand = {'conc. of A';'conc. of B'};
D(1).u = 1;
D(1).name = 'differential protein expression';
trajectories = nan(numel(t),n_data,2);
for i = 1:w*n_data
    theta_s1([4,5]) = m_gamma+randn(1,2)*sigma_gamma;
    sol = simulate_onestage(t,theta_s1,[]);
    D(1).y(1,:,i,:) = exp(sol.y+randn(numel(t),2)*sigma_noise);
    trajectories(:,i,:) = exp(sol.y);
end
for i = w*n_data+1:n_data
    theta_s2([4,5]) = m_gamma+randn(1,2)*sigma_gamma;
    sol = simulate_onestage(t,theta_s2,[]);    
    D(1).y(1,:,i,:) = exp(sol.y+randn(numel(t),2)*sigma_noise);
    trajectories(:,i,:) = exp(sol.y);
end

% test SP approximation
%test_onestage_SP(t,xi,u)
%test_onestage_SP(t,[theta_s2,2*log,u)

%%
load_plot_settings

options_plot.x_scale = 'log';
options_plot.data.plot = 'empty';
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.model.level_linewidth = 0.5;
options_plot.data.col{1} = color.data;
options_plot.model.col{1} =  'k';
options_plot.subplot_lin = 1;
options_plot.legendflag = 0;
options_plot.plainstyle = 1;
options_plot.data.kde = 1;
options_plot.marginals = 0;

for d = 1:4
options_plot.model.levelsets{1,d} = 3;
end
ODEMM_oneStage_SP
plotODEMix(D,M,xi_true,[],options_plot)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 10 4.5])
feval('print', '-dpdf','-r1000','./project/figures/forabstract.pdf');

%%
save diffgeneexp_2D_data D

D_ = D;
clear D
D(1).t = D_(1).t;
D(1).u = D_(1).u;
D(1).measurand = D_(1).measurand{1};
D(1).n_dim = 1;
D(1).name = 'differential gene expression - A';

D(1).y = D_(1).y(:,:,:,1);

D(2).t = D_(1).t;
D(2).u = D_(1).u;
D(2).measurand = D_(1).measurand{2};
D(2).n_dim = 1;
D(2).name = 'differential gene expression - B';
D(2).y = D_(1).y(:,:,:,2);

save diffgeneexp_1D_data D