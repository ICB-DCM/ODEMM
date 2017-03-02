% Visualization script for the mean pErk1/2 levels under PDL and Col I


clear all
close all
clc
load_plot_settings
load('./project/data/data_matrices_1D2D')
%%
figure('name','pErk levels - kinetic')

plot(D(1).t([1:numel(D(1).t)]), squeeze(D(1).Ey(:,:,1)), '-o',...
    'MarkerSize',2,'color',color.Lys_data); hold on;
plot(D(1).t([1:numel(D(1).t)]), squeeze(D(1).Ey(:,:,2)), '-d',...
    'MarkerSize',2,'color',color.Lys_data); hold on;
plot(D(1).t([1:numel(D(1).t)]), squeeze(D(1).Ey(:,:,3)), '-+',...
    'MarkerSize',2,'color',color.Lys_data); hold on;

plot(D(2).t([1:numel(D(1).t)]), squeeze(D(2).Ey(:,:,1)), '-o',...
    'MarkerSize',2,'color', color.ColI_data); hold on;
plot(D(2).t([1:numel(D(1).t)]), squeeze(D(2).Ey(:,:,2)), '-d',...
    'MarkerSize',2,'color', color.ColI_data); hold on;
plot(D(2).t([1:numel(D(1).t)]), squeeze(D(2).Ey(:,:,3)), '-+',...
    'MarkerSize',2,'color', color.ColI_data); hold on;
% legend('repl. 1 cond. 1','repl. 2 cond. 1','repl. 3 cond. 1',...
%     'repl. 1 cond. 2','repl. 2 cond. 2','repl. 3 cond. 2','location','southeast')
xlim([-5,max(D(1).t)+1]);
set(gca,'xtick',[0:60:120],'xticklabel',{'0','60','120'});
ylim([0,4e3])
set(gca,'TickDir','out')
xlabel('time [min]','fontsize',fs)
ylabel('pErk1/2 levels [UI]','fontsize',fs)
box off
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.1,3])
feval('print', '-dpdf','-r1000','./project/figures/pErklevels_12.pdf');
%%
figure('name','pErk levels - dose response')
plot(1:numel(D(3).u), squeeze(D(3).Ey(:,:,1)), '-o',...
    'MarkerSize',2,'color',color.Lys_data); hold on;
plot(1:numel(D(3).u), squeeze(D(3).Ey(:,:,2)), '-d',...
    'MarkerSize',2,'color',color.Lys_data); hold on;
plot(1:numel(D(3).u), squeeze(D(3).Ey(:,:,3)), '-+',...
    'MarkerSize',2,'color',color.Lys_data); hold on;
plot(1:numel(D(3).u), squeeze(D(3).Ey(:,:,4)), '-s',...
    'MarkerSize',2,'color',color.Lys_data); hold on;

plot(1:numel(D(4).u), squeeze(D(4).Ey(:,:,1)), '-o',...
    'MarkerSize',2,'color',color.ColI_data); hold on;
plot(1:numel(D(4).u), squeeze(D(4).Ey(:,:,2)), '-d',...
    'MarkerSize',2,'color',color.ColI_data); hold on;
plot(1:numel(D(4).u), squeeze(D(4).Ey(:,:,3)), '-+',...
    'MarkerSize',2,'color',color.ColI_data); hold on;
plot(1:numel(D(4).u), squeeze(D(4).Ey(:,:,4)), '-s',...
    'MarkerSize',2,'color',color.ColI_data); hold on;

xlim([1,numel(D(4).u)]);
ylim([0,2000]);

set(gca,'xtick',[1:numel(D(4).u)],'xticklabel',{'0','0.16','0.8','4','20','100','500'});
set(gca,'TickDir','out')
xlabel('NGF [ng/ml]','fontsize',fs)
ylabel('pErk1/2 levels [UI]','fontsize',fs)
box off
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.1,3])
feval('print', '-dpdf','-r1000','./project/figures/pErklevels_34.pdf');
%% average incrase
(squeeze(mean(D(4).Ey(:,:,:)))-squeeze(mean(D(3).Ey(:,:,:))))./squeeze(mean(D(3).Ey(:,:,:)))
%%
figure('name','pErk levels - dose response (TrkA)')
plot(1:numel(D(5).u), squeeze(D(5).Ey(:,:,1,2)), '-o',...
    'MarkerSize',2,'color',color.Lys_data); hold on;
plot(1:numel(D(5).u), squeeze(D(5).Ey(:,:,2,2)), '-d',...
    'MarkerSize',2,'color',color.Lys_data); hold on;
plot(1:numel(D(5).u), squeeze(D(5).Ey(:,:,3,2)), '-+',...
    'MarkerSize',2,'color',color.Lys_data); hold on;
plot(1:numel(D(5).u), squeeze(D(5).Ey(:,:,4,2)), '-s',...
    'MarkerSize',2,'color',color.Lys_data); hold on;

plot(1:numel(D(6).u), squeeze(D(6).Ey(:,:,1,2)), '-o',...
    'MarkerSize',2,'color',color.ColI_data); hold on;
plot(1:numel(D(6).u), squeeze(D(6).Ey(:,:,2,2)), '-d',...
    'MarkerSize',2,'color',color.ColI_data); hold on;
plot(1:numel(D(6).u), squeeze(D(6).Ey(:,:,3,2)), '-+',...
    'MarkerSize',2,'color',color.ColI_data); hold on;
plot(1:numel(D(6).u), squeeze(D(6).Ey(:,:,4,2)), '-s',...
    'MarkerSize',2,'color',color.ColI_data); hold on;

xlim([1,numel(D(5).u)]);
ylim([0,4000])
set(gca,'xtick',[1:numel(D(5).u)],'xticklabel',{'0','0.16','0.8','4','20','100'});
set(gca,'TickDir','out')
xlabel('NGF [ng/ml]','fontsize',fs)
ylabel('pErk1/2 levels [UI]','fontsize',fs)
box off
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.1 3])
feval('print', '-dpdf','-r1000','./project/figures/pErklevels_56.pdf');

%%
figure('name','pErk levels - dose response (Erk)')
plot(1:numel(D(7).u), squeeze(D(7).Ey(:,:,1,2)), '-o',...
    'MarkerSize',2,'color',color.Lys_data); hold on;
plot(1:numel(D(7).u), squeeze(D(7).Ey(:,:,2,2)), '-d',...
    'MarkerSize',2,'color',color.Lys_data); hold on;
plot(1:numel(D(7).u), squeeze(D(7).Ey(:,:,3,2)), '-+',...
    'MarkerSize',2,'color',color.Lys_data); hold on;
plot(1:numel(D(7).u), squeeze(D(7).Ey(:,:,4,2)), '-s',...
    'MarkerSize',2,'color',color.Lys_data); hold on;

plot(1:numel(D(8).u), squeeze(D(8).Ey(:,:,1,2)), '-o',...
    'MarkerSize',2,'color',color.ColI_data); hold on;
plot(1:numel(D(8).u), squeeze(D(8).Ey(:,:,2,2)), '-d',...
    'MarkerSize',2,'color',color.ColI_data); hold on;
plot(1:numel(D(8).u), squeeze(D(8).Ey(:,:,3,2)), '-+',...
    'MarkerSize',2,'color',color.ColI_data); hold on;
plot(1:numel(D(8).u), squeeze(D(8).Ey(:,:,4,2)), '-s',...
    'MarkerSize',2,'color',color.ColI_data); hold on;

xlim([1,numel(D(7).u)]);
set(gca,'xtick',[1:numel(D(4).u)],'xticklabel',{'0','0.16','0.8','4','20','100'});
set(gca,'TickDir','out')
xlabel('NGF [ng/ml]','fontsize',fs)
ylabel('pErk1/2 levels [UI]','fontsize',fs)
box off
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.1 3 ])
feval('print', '-dpdf','-r1000','./project/figures/pErklevels_78.pdf');

