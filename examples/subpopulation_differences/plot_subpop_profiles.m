% Visualization script for the profile likelihoods.

close all
clear all
clc

%% differences
close all
load('./results/results_subpop_TrkA.mat')
load_plot_settings

options.plot = PestoPlottingOptions();
options.plot.P.plot_type = 1;
options.plot.A.plot_type = 0;
options.plot.MS.plot_type = 0;
options.plot.P.lw = 0.8;
options.plot.boundary.mark = 0;
plotParameterProfiles(parameters,[],[],[1:16],options.plot)
legend('off')
for i = 1:16
    subplot(4,4,i)
ylabel({'likelihood'; 'ratio'},'fontsize',fs)
set(gca,'fontsize',fs)
end
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
%print('-dpdf',['./figures/subpop_profiles_1'])

plotParameterProfiles(parameters,[],[],[17:parameters.number],options.plot)
legend('off')
for i = 1:13
    subplot(4,4,i)
ylabel({'likelihood'; 'ratio'},'fontsize',fs)
set(gca,'fontsize',fs)
end
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
%print('-dpdf',['./figures/subpop_profiles_2'])


