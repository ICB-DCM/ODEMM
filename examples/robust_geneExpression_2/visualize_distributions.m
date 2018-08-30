clear all
close all
clc

load_plot_settings_robust

xgrid = linspace(1,100)';
m = 50;
C = 60;

subplot(1,4,1)
mu = m;
Sigma = C;
plot(xgrid,exp(logofmvnpdf(xgrid,mu,Sigma)),'Color',color.norm,'LineWidth',2);
xlabel('y')
ylabel('p(y)')
ylim([0,0.08])
box off
set(gca,'TickDir','out','ytick',[0,0.04,0.08])
title('normal, \varphi = (\mu,\Sigma)')

subplot(1,4,2)
delta = 0;
mu = m - sqrt(2/pi)*delta;
Sigma = C-(1-2/pi)*(delta*delta');
%plot(xgrid,exp(logofskewnormpdf(xgrid,mu,Sigma,delta)),'--','Color',color.skew_norm,'LineWidth',2);
plot(xgrid,exp(logofskewnormpdf(xgrid,mu,Sigma,delta)),'-','Color',color.skew_norm,'LineWidth',2); hold on;
delta = 12;
mu = m - sqrt(2/pi)*delta;
Sigma = C-(1-2/pi)*(delta*delta');
plot(xgrid,exp(logofskewnormpdf(xgrid,mu,Sigma,delta)),'Color',color.skew_norm_2,'LineWidth',1.5); hold on;
xlabel('y')
ylim([0,0.08])
box off
set(gca,'TickDir','out','ytick',[0,0.04,0.08])
title('skew normal')
legend('\delta = 0','\delta = 12')
legend('boxoff')

subplot(1,4,3)
mu = m;
Sigma = C;
nu = 100;
%plot(xgrid,exp(logofmvtpdf(xgrid,mu,Sigma,nu)),'--','Color',color.students_t,'LineWidth',2);
plot(xgrid,exp(logofmvtpdf(xgrid,mu,Sigma,nu)),'-','Color',color.students_t,'LineWidth',2); hold on;
nu = 2.1;
plot(xgrid,exp(logofmvtpdf(xgrid,mu,Sigma,nu)),'Color',color.students_t_2,'LineWidth',1.5); hold on;
xlabel('y')  
ylim([0,0.08])
box off
set(gca,'TickDir','out','ytick',[0,0.04,0.08])
title('Students''t')
legend('\nu = 100','\nu = 2.1')
legend('boxoff')


subplot(1,4,4)
rho = m/C;
tau = rho*m/(1-rho);
plot(xgrid,exp(logofnbinpdf(xgrid,tau,rho)),'Color',color.neg_binomial,'LineWidth',2); hold on;
xlabel('y')
ylim([0,0.08])
box off
set(gca,'TickDir','out','ytick',[0,0.04,0.08])
title('negative binomial')

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 20 5])
print('-depsc',['./figures/distributions'])

%%
figure
binwidth = 15;
k = 4;
s1 = subplot('Position',[0.1,0.25,0.4,0.7]);
h = histogram(squeeze(D(1).no_outlier(:,k,:)),'BinWidth',binwidth); hold on;
h.EdgeColor = [0.6,0.6,0.6];
h.FaceColor = color.data;
h.FaceAlpha = 1;
box off
set(gca,'TickDir','out')
ylabel('number of cells','FontSize',bigfs)
xlabel('A [counts]','FontSize',bigfs)
ylims_no = get(gca,'ylim');

s2 = subplot('Position',[0.55,0.25,0.4,0.7]);
h = histogram(squeeze(D(1).y(:,k,:)),'BinWidth',binwidth); hold on;
h.EdgeColor = [0.6,0.6,0.6];
h.FaceColor = color.data;
h.FaceAlpha = 1;
h2 = histogram(squeeze(D(1).y(:,k,D(1).ind_outlier(k,:))),'BinWidth',binwidth); hold on;
h2.EdgeColor = color.outlier;
h2.FaceColor = color.outlier;
h2.FaceAlpha = 0.4;
h2.EdgeAlpha = 0.4;
ylims = get(gca,'ylim');
ylim([0,max(ylims(2),ylims_no(2))]);
xlims = get(gca,'xlim');
box off
xlabel('A [counts]','FontSize',bigfs)
set(gca,'TickDir','out','yticklabel','')
subplot(s1);
ylim([0,max(ylims(2),ylims_no(2))]);
xlim(xlims);

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 10 4])
print('-dpdf',['./figures/outliers'])