% This script visualizes the confidence intervals for the one- and the
% two-dimensional models

close all
clear all
close all

load_plot_settings

alpha = [0.99,0.95,0.9,0.8];
delta = 0;
epsilon1 = 1e-3;
epsilon2 = 2e-3;
%% load results
load('./project/results/results_diffgeneexp_2D')
parameters2D = parameters;
clear parameters;
load('./project/results/results_diffgeneexp_1D')
parameters1D = parameters;
parameters1D_2 = parameters_2ndmode;
clear parameters parameters_2ndmode;
xi_true = [1.7,2.7,2,2,2.7,-1,-1,-1,0.5];
names = {'\lambda_{A,1}','\lambda_{A,2}','\lambda_{B,1}','\lambda_{B,2}'};
%% plot results for full distribution
figure('name','Confidence intervals');
options.PL = PestoOptions();
options.PL.mode = 'silent';
options.PL.MAP_index = 1;
lh = [];
xi_min = -3;
xi_max = 3;
ind_i = [2,3,4,5];
subplot('Position',[0.05,0.2,0.18,0.7])
for i = 1:numel(ind_i)
    for j = 1:length(alpha)
        pars_2D = getParameterConfidenceIntervals(parameters2D,alpha(j),options.PL);
        pars_2D.CI.PL = max(min(pars_2D.CI.PL,xi_max),xi_min);
        fill(pars_2D.CI.PL(ind_i(i),[1,2,2,1]),(i*(-2-delta))+(1.5+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_full(j,:)); hold on;
        %fill((pars_2D.MS.par(ind_i(i))+epsilon1*[-1,+1,+1,-1]),(i*(-2-delta))+1*[+1,+1,-1,-1],'b','facecolor',color.CI_full(end,:)); hold on;
        if i == numel(ind_i)
            plot([xi_true(ind_i(i)), xi_true(ind_i(i))],[((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2+0.5,(((i-1)*(-2-delta)-1)+((i)*(-2-delta)+0))/2+0.5],'-','color',color.true,'LineWidth',0.8);
        else
            plot([xi_true(ind_i(i)), xi_true(ind_i(i))],[((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2+0.5,(((i-1)*(-2-delta)-1)+((i)*(-2-delta)+0))/2+0.5],'-','color',color.true,'LineWidth',0.8);
        end
    end
    plot([pars_2D.MS.par(ind_i(i)),pars_2D.MS.par(ind_i(i))],(i*(-2-delta))+[+1,-1],'LineWidth',1.5,'color',color.CI_full(end,:)); hold on;
    if i < numel(ind_i)
        plot([-6,6],((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2*[1,1]+0.5,'-','color',0.8*[1,1,1]);
    end
end
set(gca,'xscale','lin',...
    'xlim',[1.95,2.052],...
    'xtick',[1.95,2,2.05,2.65,2.7,2.75],...
    'ytick',[],...
    'yticklabel',{},...
    'ylim',[-9,-1],...
    'tickdir','out',...
    'fontsize',6);
box off

subplot('Position',[0.27,0.2,0.18,0.7])
for i = 1:numel(ind_i)
    for j = 1:length(alpha)
        pars_2D = getParameterConfidenceIntervals(parameters2D,alpha(j),options.PL);
        pars_2D.CI.PL = max(min(pars_2D.CI.PL,xi_max),xi_min);
        fill(pars_2D.CI.PL(ind_i(i),[1,2,2,1]),(i*(-2-delta))+(1.5+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_full(j,:)); hold on;
        %fill((pars_2D.MS.par(ind_i(i))+epsilon1*[-1,+1,+1,-1]),(i*(-2-delta))+1*[+1,+1,-1,-1],'b','facecolor',color.CI_full(end,:)); hold on;
    end
    plot([pars_2D.MS.par(ind_i(i)),pars_2D.MS.par(ind_i(i))],(i*(-2-delta))+[+1,-1],'LineWidth',1.5,'color',color.CI_full(end,:)); hold on;
    if i == numel(ind_i)
        plot([xi_true(ind_i(i)), xi_true(ind_i(i))],[((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2+0.5,(((i-1)*(-2-delta)-1)+((i)*(-2-delta)+0))/2+0.5],'-','color',color.true,'LineWidth',0.8);
    else
        plot([xi_true(ind_i(i)), xi_true(ind_i(i))],[((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2+0.5,(((i-1)*(-2-delta)-1)+((i)*(-2-delta)+0))/2+0.5],'-','color',color.true,'LineWidth',0.8);
    end
    if i < numel(ind_i)
        plot([-6,6],((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2*[1,1]+0.5,'-','color',0.8*[1,1,1]);
    end
end
set(gca,'xscale','lin',...
    'xlim',[2.64,2.75],...
    'xtick',[1.95,2,2.05,2.65,2.7,2.75],...
    'ytick',[],...
    'ylim',[-9,-1],...
    'tickdir','out',...
    'yticklabel',{},'fontsize',6);
box off

%% plot results for marginal distributions
subplot('Position',[0.55,0.2,0.18,0.7])
for i = 1:numel(ind_i)
    for j = 1:length(alpha)
        options.PL.MAP_index = 1;
        pars_1D = getParameterConfidenceIntervals(parameters1D,alpha(j),options.PL);
        pars_1D.CI.PL = max(min(pars_1D.CI.PL,xi_max),xi_min);
        options.PL.MAP_index = 43;
        pars_1D_2 = getParameterConfidenceIntervals(parameters1D_2,alpha(j),options.PL);
        if ind_i(i) == 4 || ind_i(i) == 5
            fill(pars_1D_2.CI.PL(ind_i(i),[1,2,2,1]),(i*(-2-delta))+(1.5+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(j,:)); hold on;
        end
        fill(pars_1D.CI.PL(ind_i(i),[1,2,2,1]),(i*(-2-delta))+(1.5+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(j,:)); hold on;
        %fill((pars_1D.MS.par(ind_i(i),3)+epsilon1*[-1,+1,+1,-1]),(i*(-2-delta))+1*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(end,:)); hold on;
        
    end
    plot([pars_1D.MS.par(ind_i(i)),pars_1D.MS.par(ind_i(i))],(i*(-2-delta))+[+1,-1],'LineWidth',1.5,'color',color.CI_marg(end,:)); hold on;
    if i == numel(ind_i)
        plot([xi_true(ind_i(i)), xi_true(ind_i(i))],[((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2+0.5,(((i-1)*(-2-delta)-1)+((i)*(-2-delta)+0))/2+0.5],'-','color',color.true,'LineWidth',0.8);
    else
        plot([xi_true(ind_i(i)), xi_true(ind_i(i))],[((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2+0.5,(((i-1)*(-2-delta)-1)+((i)*(-2-delta)+0))/2+0.5],'-','color',color.true,'LineWidth',0.8);
    end
    if i < numel(ind_i)
        plot([-6,6],((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2*[1,1]+0.5,'-','color',0.8*[1,1,1]);
    end
end
set(gca,'xscale','lin',...
    'xlim',[1.95,2.052],...
    'xtick',[1.95,2,2.05,2.65,2.7,2.75],...
    'ytick',[],...
    'yticklabel',{},...
    'ylim',[-9,-1],...
    'tickdir','out',...
    'fontsize',6);
box off
subplot('Position',[0.77,0.2,0.18,0.7])
for i = 1:numel(ind_i)
    for j = 1:length(alpha)
        options.PL.MAP_index = 1;
        pars_1D = getParameterConfidenceIntervals(parameters1D,alpha(j),options.PL);
        pars_1D.CI.PL = max(min(pars_1D.CI.PL,xi_max),xi_min);
        options.PL.MAP_index = 43;
        pars_1D_2 = getParameterConfidenceIntervals(parameters1D_2,alpha(j),options.PL);
        if ind_i(i) == 4 || ind_i(i) == 5
            fill(pars_1D_2.CI.PL(ind_i(i),[1,2,2,1]),(i*(-2-delta))+(1.5+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(j,:)); hold on;
        end
        fill(pars_1D.CI.PL(ind_i(i),[1,2,2,1]),(i*(-2-delta))+(1.5+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(j,:)); hold on;
        %fill((pars_1D.MS.par(ind_i(i),3)+epsilon1*[-1,+1,+1,-1]),(i*(-2-delta))+1*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(end,:)); hold on;
        if i == numel(ind_i)
            plot([xi_true(ind_i(i)), xi_true(ind_i(i))],[((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2+0.5,(((i-1)*(-2-delta)-1)+((i)*(-2-delta)+0))/2+0.5],'-','color',color.true,'LineWidth',0.8);
        else
            plot([xi_true(ind_i(i)), xi_true(ind_i(i))],[((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2+0.5,(((i-1)*(-2-delta)-1)+((i)*(-2-delta)+0))/2+0.5],'-','color',color.true,'LineWidth',0.8);
        end
    end
    plot([pars_1D.MS.par(ind_i(i)),pars_1D.MS.par(ind_i(i))],(i*(-2-delta))+[+1,-1],'LineWidth',1.5,'color',color.CI_marg(end,:)); hold on;
    
    if i < numel(ind_i)
        plot([-6,6],((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2*[1,1]+0.5,'-','color',0.8*[1,1,1]);
    end
end
set(gca,'xscale','lin',...
    'xlim',[2.64,2.75],...
    'xtick',[1.95,2,2.05,2.65,2.7,2.75],...
    'ytick',[],...
    'ylim',[-9,-1],...
    'tickdir','out',...
    'yticklabel',{},'fontsize',fs);
box off
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 8.35 2.2])
print('-depsc','./project/figures/CI_manuscript');

% %% Uncertainty Figure Supplement
% figure('name','Confidence Interval Supplement');
% xi_min = -1.5;
% xi_max = -0.7;
% ind = [6,7,8];
% subplot('Position',[0.1,0.2,0.5,0.75]);
% for i = 1:numel(ind)
%     for j = 1:length(alpha)
%         pars_2D = getParameterConfidenceIntervals(parameters2D,alpha(j));
%         options.PL.MAP_index = 1;
%         pars_1D = getParameterConfidenceIntervals(parameters1D,alpha(j),options.PL);
%         pars_2D.CI.PL = max(min(pars_2D.CI.PL,xi_max),xi_min);
%         pars_1D.CI.PL = max(min(pars_1D.CI.PL,xi_max),xi_min);
%         fill(pars_2D.CI.PL(ind(i),[1,2,2,1]),(i*(-2-delta))+(0.45+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_full(j,:)); hold on;
%         fill(pars_1D.CI.PL(ind(i),[1,2,2,1]),(i*(-2-delta)-1)+(0.45+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(j,:)); hold on;
%         %fill((pars_2D.MS.par(ind(i))+epsilon1*[-1,+1,+1,-1]),(i*(-2-delta))+0.45*[+1,+1,-1,-1],'b','facecolor',color.CI_full(end,:)); hold on;
%         %fill((pars_1D.MS.par(ind(i))+epsilon1*[-1,+1,+1,-1]),(i*(-2-delta)-1)+0.45*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(end,:)); hold on;
%     end
%     plot([pars_1D.MS.par(ind(i)),pars_1D.MS.par(ind(i))],(i*(-2-delta)-1)+0.45*[+1,-1],'LineWidth',1.5,'color',color.CI_marg(end,:)); hold on;
%     plot([pars_2D.MS.par(ind(i)),pars_2D.MS.par(ind(i))],(i*(-2-delta))+0.45*[+1,-1],'LineWidth',1.5,'color',color.CI_full(end,:)); hold on;
%     plot([xi_true(ind(i)), xi_true(ind(i))],[((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2,(((i-1)*(-2-delta)-1)+((i)*(-2-delta)+0))/2],'-','color',color.true,'LineWidth',0.8);
%     
%     if i < numel(ind)
%         plot([-2,0.6],((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2*[1,1],'-','color',0.8*[1,1,1]);
%     end
% end
% set(gca,'xscale','lin',...
%     'xlim',[-1.3,-0.7],...
%     'xtick',[-1.3:.1:-0.7],...
%     'ytick',[],...
%     'tickdir','out',...
%     'ylim',[-7.5,-1.5],...
%     'yticklabel',{},'fontsize',fs);
% %set(gca,'Position',[0.2,0.15,0.75,0.83])
% box off
% %set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 6 4])
% %print('-depsc','./project/figures/CI_rest');
% %%
% subplot('Position',[0.7,0.2,0.2,0.25]);
% xi_min = 0;
% xi_max = 1;
% ind=9; %weight
% i=1;
% for j = 1:length(alpha)
%     pars_2D = getParameterConfidenceIntervals(parameters2D,alpha(j));
%     pars_1D = getParameterConfidenceIntervals(parameters1D,alpha(j),options.PL);
%     pars_2D.CI.PL = max(min(pars_2D.CI.PL,xi_max),xi_min);
%     pars_1D.CI.PL = max(min(pars_1D.CI.PL,xi_max),xi_min);
%     fill(pars_2D.CI.PL(ind(i),[1,2,2,1]),(i*(-2-delta))+(0.45+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_full(j,:)); hold on;
%     fill(pars_1D.CI.PL(ind(i),[1,2,2,1]),(i*(-2-delta)-1)+(0.45+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(j,:)); hold on;
%     %fill((pars_2D.MS.par(ind(i))+epsilon1*[-1,+1,+1,-1]),(i*(-2-delta))+0.45*[+1,+1,-1,-1],'b','facecolor',color.CI_full(end,:)); hold on;
%     %fill((pars_1D.MS.par(ind(i))+epsilon1*[-1,+1,+1,-1]),(i*(-2-delta)-1)+0.45*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(end,:)); hold on;
% end
% plot([pars_1D.MS.par(ind(i)),pars_1D.MS.par(ind(i))],(i*(-2-delta)-1)+0.45*[+1,-1],'LineWidth',1.5,'color',color.CI_marg(end,:)); hold on;
% plot([pars_2D.MS.par(ind(i)),pars_2D.MS.par(ind(i))],(i*(-2-delta))+0.45*[+1,-1],'LineWidth',1.5,'color',color.CI_full(end,:)); hold on;
% plot([xi_true(ind(i)), xi_true(ind(i))],[((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2,(((i-1)*(-2-delta)-1)+((i)*(-2-delta)+0))/2],'-','color',color.true,'LineWidth',0.8);
% 
% set(gca,'xscale','lin',...
%     'xlim',[0.4,0.6],...
%     'xtick',[0.4,0.5,0.6],...
%     'ytick',[],...
%     'tickdir','out',...
%     'ylim',[-3.5,-1.5],...
%     'yticklabel',{},...
%     'fontsize',6);
% box off
% xi_min = 1.2;
% xi_max = 2;
% subplot('Position',[0.7,0.7,0.2,0.25]);
% ind=1;
% for j = 1:length(alpha)
%     pars_2D = getParameterConfidenceIntervals(parameters2D,alpha(j));
%     options.PL.MAP_index = 1;
%     pars_1D = getParameterConfidenceIntervals(parameters1D,alpha(j),options.PL);
%     pars_2D.CI.PL = max(min(pars_2D.CI.PL,xi_max),xi_min);
%     pars_1D.CI.PL = max(min(pars_1D.CI.PL,xi_max),xi_min);
%     fill(pars_2D.CI.PL(ind(i),[1,2,2,1]),(i*(-2-delta))+(0.45+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_full(j,:)); hold on;
%     fill(pars_1D.CI.PL(ind(i),[1,2,2,1]),(i*(-2-delta)-1)+(0.45+j*0.05)*0.5*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(j,:)); hold on;
%     %fill((pars_2D.MS.par(ind(i))+epsilon1*[-1,+1,+1,-1]),(i*(-2-delta))+0.45*[+1,+1,-1,-1],'b','facecolor',color.CI_full(end,:)); hold on;
%     %fill((pars_1D.MS.par(ind(i))+epsilon1*[-1,+1,+1,-1]),(i*(-2-delta)-1)+0.45*[+1,+1,-1,-1],'b','facecolor',color.CI_marg(end,:)); hold on;
% end
% plot([pars_1D.MS.par(ind(i)),pars_1D.MS.par(ind(i))],(i*(-2-delta)-1)+0.45*[+1,-1],'LineWidth',1.5,'color',color.CI_marg(end,:)); hold on;
% plot([pars_2D.MS.par(ind(i)),pars_2D.MS.par(ind(i))],(i*(-2-delta))+0.45*[+1,-1],'LineWidth',1.5,'color',color.CI_full(end,:)); hold on;
% plot([xi_true(ind(i)), xi_true(ind(i))],[((i*(-2-delta)-1)+((i+1)*(-2-delta)+0))/2,(((i-1)*(-2-delta)-1)+((i)*(-2-delta)+0))/2],'-','color',color.true,'LineWidth',0.8);
% 
% set(gca,'xscale','lin',...
%     'xlim',[1.4,2],...
%     ...'xtick',[0.4,0.5,0.6],...
%     'ytick',[],...
%     'tickdir','out',...
%     'ylim',[-3.5,-1.5],...
%     'yticklabel',{},...
%     'fontsize',6);
% box off
% 
% set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 13 5])
% print('-depsc','./project/figures/CI_supplement');
% 
% %% check convergence
% conv(1) = sum(2*(parameters2D.MS.logPost-parameters2D.MS.logPost(1))>-icdf('chi2',0.95,1));
% conv(2) = sum(2*(parameters1D.MS.logPost-parameters1D.MS.logPost(1))>-icdf('chi2',0.95,1));



