% This script generates the artificial data of a conversion process

close all
clear all
clc

n_data = 1e3;
%% Data generation
t = [0,0.1,0.5,1,2]; %hours
xi_true = [-0.1,0.1,-0.45,-0.2,-1,-1.8,0.7];

theta_s1 = log(10)*xi_true([1,3,4]); % log
theta_s2 = log(10)*xi_true([2,3,4]); % log 
sigma_k3 = 10.^(xi_true(5)); % lin

w = xi_true(end); % lin
sigma_noise = 10.^xi_true(end-1); % lin

m_k3 = theta_s1(3); 

D(1).y = nan(1,numel(t),n_data,1);
D(1).t = t;
D(1).n_dim = 1;
D(1).measurand = 'conc. of B';
D(1).u = 1;
D(1).name = 'conversion reaction';
trajectories = nan(numel(t),n_data);
for i = 1:w*n_data
    i
    theta_s1(3) = m_k3+randn(1)*sigma_k3;
    sol = simulate_CR_log(t,theta_s1,[]);
    D(1).y(1,:,i,:) = exp(sol.y'+randn(size(sol.y'))*sigma_noise);
    trajectories(:,i) = exp(sol.y);
    
end
for i = w*n_data+1:n_data
    i
    theta_s2(3) = (m_k3+randn(1)*sigma_k3);
    sol = simulate_CR_log(t,theta_s2,[]);
    D(1).y(1,:,i,:) = exp(sol.y'+randn(size(sol.y'))*sigma_noise);
    trajectories(:,i) = exp(sol.y);
end
options_plot.data.bins = 20;
plotODEMix(D)

%% check Sigma Point approximation
test_SP_CR(t,[theta_s1,2*log(sigma_k3)],[])
test_SP_CR(t,[theta_s2,2*log(sigma_k3)],[])
%% save data
save project/data/conversionprocess_data D xi_true

%save(['conversionprocess_data_' num2str(n_data)],'D','xi_true')
