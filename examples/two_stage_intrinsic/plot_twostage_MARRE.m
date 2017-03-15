% Visualization script for the two stage gene expression example with
% intrinsic noise.

clear all
close all
clc

load_plot_settings
%% Fit for Moment Approximation
options_plot.x_scale = 'lin';
options_plot.data.plot = 'filled';
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.data.fill_col{1} = color.data;
options_plot.model.col{1} = color.SP;
options_plot.boundaries(1).y_min = 0;
options_plot.boundaries(1).y_max = 2;
options_plot.model.lw = 1;
options_plot.model.ls = '-';

load('./project/results/geneExp_MA','M','D','parameters','conditions')
xi = parameters.MS.par(:,1);
fh = plotODEMix(D,M,xi,[],options_plot);

%% Evaluation of BIC, CPU time and convergence for Moment Approximation
parameters.MS.BIC = -2*parameters.MS.logPost+ log(sum(sum(~isnan(D(1).y))))*parameters.number;
BICs(1) = parameters.MS.BIC(1);
np(1) = parameters.number;
t_cpu(1) = nanmean(parameters.MS.t_cpu(parameters.MS.t_cpu>0));
t_cpus(1,1:numel(parameters.MS.t_cpu(parameters.MS.t_cpu>0))) = parameters.MS.t_cpu(parameters.MS.t_cpu>0);
conv(1) = sum(2*(parameters.MS.logPost-parameters.MS.logPost(1))>-icdf('chi2',0.95,1))/parameters.MS.n_starts;
p1 = parameters;

%% Evaluation of BIC, CPU time and convergence for RRE and fit
load('./project/results/geneExp_RRE','M','D','parameters','conditions')
xi = parameters.MS.par(:,1);
parameters.MS.BIC = -2*parameters.MS.logPost+ log(sum(sum(~isnan(D(1).y))))*parameters.number;
BICs(2) = parameters.MS.BIC(1);
np(2) = parameters.number;
t_cpu(2) = nanmean(parameters.MS.t_cpu(parameters.MS.t_cpu>0));
t_cpus(2,1:numel(parameters.MS.t_cpu(parameters.MS.t_cpu>0))) = parameters.MS.t_cpu(parameters.MS.t_cpu>0);
t_cpus = t_cpus';
conv(2) = sum(2*(parameters.MS.logPost-parameters.MS.logPost(1))>-icdf('chi2',0.95,1))/parameters.MS.n_starts;
p2 = parameters;
load('./project/results/geneExp_RRE_2ndMode','parameters_2')
p2_2 = parameters_2;

options_plot.model.col{1} = color.RRE;
options_plot.hold_on = 1;
options_plot.model.ls = '--';
D_ = D;
D_.y = [];
plotODEMix(D_,M,xi,[],options_plot,[],fh);
ah=axes('position',[.13,.11,0.7813,0],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.6,'Color','k');
ah=axes('position',[.13,.11,0,0.85],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.6,'Color','k');
s = subplot(1,numel(D(1).t),1);
set(gca,'xtick',[0,1,2]);
set(gca,'tickdir','out')
for i = 1:numel(D(1).t)
    s = subplot(1,numel(D(1).t),i);
    ylim([0,0.4]);
end
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7 3])
set(gcf,'renderer','opengl');
feval('print', '-dpdf','-r1000',['./project/figures/datafit_MARRE']);

%% Convergence plot
figure('name','Convergence')
plot(1:p1.MS.n_starts,p1.MS.logPost,'o','MarkerSize',3,'color',color.SP); hold on;
plot(1:p2.MS.n_starts,p2.MS.logPost,'o','MarkerSize',3,'color',color.RRE);
ylabel('log-likelihood','fontsize',fs)
xlabel('sorted optimizer start','fontsize',fs)
box off
set(gca,'Tickdir','out')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 6 4])
feval('print', '-dpdf',['./project/figures/conv']);


figure('name','Convergence Zoom')
plot(1:p1.MS.n_starts,p1.MS.logPost,'o','MarkerSize',3,'color',color.SP); hold on;
plot(1:p2.MS.n_starts,p2.MS.logPost,'o','MarkerSize',3,'color',color.RRE)
ylim([3610,3630])
xlim([0,61])
box off
set(gca,'ytick',[3610,3630],'xtick',[0,20,40,60]),
ylabel('')
xlabel('sorted optimizer start','fontsize',fs)
set(gca,'Tickdir','out')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4 2])
feval('print', '-dpdf',['./project/figures/conv_zoom']);

legend('MA','RRE','Location','SouthWest')
legend('boxoff')
set(gca,'ytick',[],'xtick',[]);
xlabel('');
ylim([0,1])
xlim([0,1])
axis off
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4 2])
feval('print', '-dpdf',['./project/figures/conv_leged']);

conv(1) = sum(2*(p1.MS.logPost-p1.MS.logPost(1))>-icdf('chi2',0.95,1));
conv(2) = sum(2*(p2.MS.logPost-p2.MS.logPost(1))>-icdf('chi2',0.95,1));

disp(['convergence for MA: ' num2str(conv(1))]);
disp(['convergence for RRE: ' num2str(conv(2))]);
%% cpu time
figure('name','CPU time')
t_cpus(t_cpus==0) = NaN;
boxplot(t_cpus(:,[2,1]),{'RRE','MA'},'Color','k','symbol','o',...
    'outliersize',1)

h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),color.SP,'FaceAlpha',.7);
patch(get(h(2),'XData'),get(h(2),'YData'),color.RRE,'FaceAlpha',.7);
set(gca,'tickdir','out','ylim',[0,75]);
hold on;
xdata1 = get(h(2),'XData');
xdata2 = get(h(1),'XData');
box off
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4 4])
ylabel('cpu time [sec]','fontsize',fs)
feval('print', '-dpdf',['./project/figures/tcpus']);

%% visualize BIC
figure('name','BIC')
plot(1,BICs(2)-min(BICs),'.','MarkerSize',10,'Color',color.RRE); hold on;
plot(2,BICs(1)-min(BICs),'.','MarkerSize',10,'Color',color.SP)
box off
ylabel('BIC - min(BIC)','fontsize',fs)
xlim([0.9,2.1])
ylim([-5,80])
set(gca,'xtick',[1,2],'xticklabel',{'RRE','MA'},'tickdir','out');
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2.5 2.5])
feval('print', '-dpdf',['./project/figures/BICs']);

disp(['BIC-min(BIC) for MA: ' num2str(BICs(1)-min(BICs))]);
disp(['BIC-min(BIC) for RRE: ' num2str(BICs(2)-min(BICs))]);

%% Profiles
xi_true = log10([10,10,20,1,5,0.1]);
xi_true(7) = 0.5;
names = {'k_1','k_{2,1}','k_{2,2}','k_3','k_4','k_5','w_1'};
xlimits{1} = [8,12];
xlimits{2} = [8,12];
xlimits{3} = [15,25];
xlimits{4} = [0.09,2];
xlimits{5} = [2,10];
xlimits{6} = [0.09,2];
xlimits{7} = [0.4,0.6];
xticklabels{1} = {'8','10','12'};
xticklabels{2} = {'8','10','12'};
xticklabels{3} = {'15','20','25'};
xticklabels{4} = {'0.09','1','2'};
xticklabels{5} = {'2','6','10'};
xticklabels{6} = {'0.09','1','2'};
xticklabels{7} = {'0.4','0.5','0.6'};
% %
% figure('name','Profiles MA vs RRE')
% pc = 1;
% for s = 1:7
%     subplot(2,4,pc);
%     pc = pc+1;
%     if pc == 4
%         pc = pc+1;
%     end
%     if s < 7
%         plot(10.^p1.P(s).par(s,:),p1.P(s).R,'Color',color.SP,'LineWidth',1.3); hold on;
%         plot(10.^p2.P(s).par(s,:),p2.P(s).R,'Color',color.RRE,'LineWidth',1.3); hold on;
%         if s == 4 || s == 6
%             plot(10.^p2_2.P(s).par(s,:),p2_2.P(s).R,'Color',color.RRE,'LineWidth',1.3); hold on;
%         end
%         plot([10.^xi_true(s),10.^xi_true(s)],[0,1],'Color',color.true,'LineWidth',1.3); hold on;
%     else
%         plot(p1.P(s).par(s,:),p1.P(s).R,'Color',color.SP,'LineWidth',1.3); hold on;
%         plot([xi_true(s),xi_true(s)],[0,1],'Color',color.true,'LineWidth',1.3); hold on;
%     end
%     ylim([0,1.01])
%     xlim(xlimits{s});
%     ylabel('likelihood ratio','fontsize',fs)
%     %set(gca,'xtick',log10(linspace(10^xlimits{s}(1),10^xlimits{s}(2),3)),'xticklabel',xticklabels{s},...
%     %    'tickdir','out')
%     set(gca,'xtick',linspace(xlimits{s}(1),xlimits{s}(2),3),'xticklabel',xticklabels{s},...
%         'tickdir','out')
%     box off
%     xlabel(names{s},'fontsize',fs)
% end
%%
figure('name','Profiles MA vs RRE 1')
pc = 1;
positions = [0.1 0.25  0.15 0.7;...
    0.34 0.25 0.15 0.7;...
    0.58 0.25 0.15 0.7;...
    0.82 0.25 0.15 0.7]
for s = [1,2,3,5]
    subplot('Position',positions(pc,:));
    pc = pc+1;
    if s < 7
        if s==2
            plot(10.^p1.P(3).par(3,:),p1.P(s).R,'Color',color.SP,'LineWidth',1.3); hold on;
        elseif s==3
            plot(10.^p1.P(2).par(2,:),p1.P(s).R,'Color',color.SP,'LineWidth',1.3); hold on
        else
            plot(10.^p1.P(s).par(s,:),p1.P(s).R,'Color',color.SP,'LineWidth',1.3); hold on
        end
        plot(10.^p2.P(s).par(s,:),p2.P(s).R,'--','Color',color.RRE,'LineWidth',1.3); hold on;
        if s == 4 || s == 6
            plot(10.^p2_2.P(s).par(s,:),p2_2.P(s).R,'--','Color',color.RRE,'LineWidth',1.3); hold on;
        end
        plot([10.^xi_true(s),10.^xi_true(s)],[0,1],'Color',color.true,'LineWidth',1.3); hold on;
    else
        plot(p1.P(s).par(s,:),p1.P(s).R,'Color',color.SP,'LineWidth',1.3); hold on;
        plot([xi_true(s),xi_true(s)],[0,1],'Color',color.true,'LineWidth',1.3); hold on;
    end
    ylim([0,1.01])
    xlim(xlimits{s});
    ylabel('likelihood ratio','fontsize',fs)
    set(gca,'xtick',linspace(xlimits{s}(1),xlimits{s}(2),3),'xticklabel',xticklabels{s},...
        'tickdir','out')
    box off
    xlabel(names{s},'fontsize',fs)
end

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 12 3])
feval('print', '-dpdf',['./project/figures/profiles_1']);

%%
figure('name','Profiles MA vs RRE 2')
pc = 1;
positions = [0.1 0.25  0.15 0.7;...
    0.26 0.25 0.15 0.7;...
    0.46 0.25 0.15 0.7;...
    0.62 0.25 0.15 0.7
    0.82 0.25 0.15 0.7];

for s = [4,6]
    for i = 1:2
        subplot('Position',positions(pc,:));
        pc = pc+1;
        plot(10.^p1.P(s).par(s,:),p1.P(s).R,'Color',color.SP,'LineWidth',1.3); hold on;
        plot(10.^p2.P(s).par(s,:),p2.P(s).R,'--','Color',color.RRE,'LineWidth',1.3); hold on;
        plot(10.^p2_2.P(s).par(s,:),p2_2.P(s).R,'--','Color',color.RRE,'LineWidth',1.3); hold on;
        plot([10.^xi_true(s),10.^xi_true(s)],[0,1],'Color',color.true,'LineWidth',1.3); hold on;
        ylim([0,1.01])
        xlim(xlimits{s});
        ylabel('likelihood ratio','fontsize',fs)
        set(gca,'xtick',[0.09,0.1,0.11],...
            'tickdir','out')
        box off
        xlim([0.09,0.112]);
        xlabel(names{s},'fontsize',fs)
        if i == 2
            set(gca,'ytick',[],'yticklabel','',...
                'tickdir','out','xtick',[0.5,1,2])
            xlim([0.41,2]);
            ylabel('');
        end
    end
end
s=7;
subplot('Position',positions(5,:));
plot(1-p1.P(s).par(s,:),p1.P(s).R,'Color',color.SP,'LineWidth',1.3); hold on;
plot(p2.P(s).par(s,:),p2.P(s).R,'--','Color',color.RRE,'LineWidth',1.3); hold on;
plot([xi_true(s),xi_true(s)],[0,1],'Color',color.true,'LineWidth',1.3); hold on;
ylim([0,1.01])
xlim(xlimits{s});
ylabel('likelihood ratio','fontsize',fs)
set(gca,'xtick',linspace(xlimits{s}(1),xlimits{s}(2),3),'xticklabel',xticklabels{s},...
    'tickdir','out')
box off
xlabel(names{s},'fontsize',fs)


set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 15 3])
feval('print', '-dpdf',['./project/figures/profiles_2']);




