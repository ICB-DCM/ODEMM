% Visualization script for the profile likelihoods for the final model 

close all
clear all
clc

%% differences
load('./project/results/results_ECM_finalmodel.mat')
close all
load_plot_settings 


parameters.name = {'log_{10}(k_1)','log_{10}(k_2)','log_{10}(k_4)','log_{10}(k_{5})',...% parameter for simulation
        'log_{10}(\beta_{k_3[TrkA]_{0,1}})','log_{10}(\beta_{k_3TrkA_{0,2}})',...
        'log_{10}(\beta_{c_P^{1,2}[Erk]_{0,1}})',...
        'log_{10}(\sigma_{T,1})','log_{10}(\sigma_{T,2})','log_{10}(\sigma_{E})','log_{10}(\sigma_{TE})'...
        'log_{10}(c_P^{1,2})','log_{10}(c_T^{5,6}/k_3)','log_{10}(o_T^{5,6})',... % scaling and offset
        'log_{10}(c_E^{7,8}/c_P^{5,6})','log_{10}(o_E^{7,8})',...
        'log_{10}(c_P^{3,4})','log_{10}(o_P^{3,4})',...
        'log_{10}(c_P^{5,6})','log_{10}(o_P^{5,6})',...
        'log_{10}(c_P^{7,8})','log_{10}(o_P^{7,8})',...
        'log_{10}(\sigma_{P,noise}^{1,2})','log_{10}(\sigma_{P,noise}^{3,4})',... %measurement noise
        'log_{10}(\sigma_{T,noise}^{5,6})','log_{10}(\sigma_{P,noise}^{5,6})',...
        'log_{10}(\sigma_{E,noise}^{7,8})','log_{10}(\sigma_{P,noise}^{7,8})',...
        'log_{10}(w)','log_{10}(\kappa_{k_5})','log_{10}(\kappa_{k_3[TrkA]_0})','log_{10}(\kappa_{c_{P}^{1,2}Erk})'};

options.plot.A.plot_type = 1;
options.plot.A.lw = 1;
options.plot.MS.lw = 0.8;
options.plot.MS.plot_type = 1;
options.plot.P.lw = 0.8;
options.plot.A.lw = 0.8;


plotParameterProfiles(parameters,[],[],[1:16],options.plot)
subplot(6,6,5)
xlim([-4,-3.5]);
subplot(6,6,6)
xlim([-1.5,-0.5]);
subplot(6,6,8)
xlim([0.6,0.9]);
subplot(6,6,9)
xlim([-0.6,-0.45]);
subplot(6,6,29)
xlim([-0.55,-0.51]);

for i = 32:-1:1
    subplot(6,6,i)
ylabel({'likelihood'; 'ratio'},'fontsize',fs)
set(gca,'fontsize',fs)
end

 legend('maximum likelihood estimate','likelihood ratio derived from profile likelihood','likelihood ratio derived from local approximation at maximum likelihood estimate')
  legend('boxoff','fontsize',fs)

 %%
 set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 14])
 print('-dpdf',['./project/figures/ECM_profiles'])

%%
plotParameterProfiles(parameters,[],[],[1:16],options.plot)
subplot(4,4,5)
xlim([-4,-3.5]);
subplot(4,4,6)
xlim([-1.5,-0.5]);
subplot(4,4,8)
xlim([0.6,0.9]);
subplot(4,4,9)
xlim([-0.6,-0.45]);

legend('off')
for i = 1:16
    subplot(4,4,i)
ylabel({'likelihood'; 'ratio'},'fontsize',fs)
set(gca,'fontsize',fs)
end
 %%
 set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 ])
 print('-dpdf',['./project/figures/ECM_profiles_1'])
%%
plotParameterProfiles(parameters,[],[],[17:parameters.number],options.plot)

subplot(4,4,13)
xlim([-0.55,-0.51]);

legend('off')
for i = 1:16
    subplot(4,4,i)
ylabel({'likelihood'; 'ratio'},'fontsize',fs)
set(gca,'fontsize',fs)
end
%%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
 print('-dpdf',['./project/figures/ECM_profiles_2'])
%%