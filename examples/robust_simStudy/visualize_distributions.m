clear all
close all
clc

load_plot_settings_robust

fh = figure;
xgrid = linspace(1,100)';
m = 50;
C = 60;

s1 = subplot(1,4,1);
mu = m;
Sigma = C;
plot(xgrid,exp(logofmvnpdf(xgrid,mu,Sigma)),'Color',...
    color.norm,'LineWidth',2);
xlabel('$\bar{y}$','interpreter','latex')
ylabel('probability density/mass')
ylim([0,0.08])
box off
s1.FontSize = fs;
s1.TickDir='out';   
s1.Position=[0.08,0.2,0.18,0.8];

s2 = subplot(1,4,2);
delta = 0;
mu = m - sqrt(2/pi)*delta;
Sigma = C-(1-2/pi)*(delta*delta');
plot(xgrid,exp(logofskewnormpdf(xgrid,mu,Sigma,delta)),'-',...
    'Color',color.skew_norm,'LineWidth',2); hold on;
delta = 12;
mu = m - sqrt(2/pi)*delta;
Sigma = C-(1-2/pi)*(delta*delta');
plot(xgrid,exp(logofskewnormpdf(xgrid,mu,Sigma,delta)),...
    'Color',color.skew_norm_2,'LineWidth',1.5); hold on;
xlabel('$\bar{y}$','interpreter','latex')
ylim([0,0.08])
box off
s2.FontSize = fs;
s2.TickDir='out'; 
s2.YTickLabel = '';
s2.Position =  [0.3,0.2,0.18,0.8];

hl1=legend('\delta = 0','\delta = 12');
legend('boxoff')
hl1.Position=[0.32,0.9,0.16,0.05];

s3 = subplot(1,4,3);
mu = m;
Sigma = C;
nu = 100;
plot(xgrid,exp(logofmvtpdf(xgrid,mu,Sigma,nu)),'-',...
    'Color',color.students_t,'LineWidth',2); hold on;
nu = 2.1;
plot(xgrid,exp(logofmvtpdf(xgrid,mu,Sigma,nu)),'Color',...
    color.students_t_2,'LineWidth',1.5); hold on;
xlabel('$\bar{y}$','interpreter','latex')
ylim([0,0.08])
box off
s3.FontSize = fs;
s3.TickDir='out'; 
s3.YTickLabel = '';
s3.Position = [0.52,0.2,0.18,0.8];
hl2 = legend('\nu = 100','\nu = 2.1');
hl2.Position=[0.54,0.9,0.16,0.05];
legend('boxoff')

s4=subplot(1,4,4);
rho = m/C;
tau = rho*m/(1-rho);
plot(xgrid,exp(logofnbinpdf(xgrid,tau,rho)),'Color',....
    color.neg_binomial,'LineWidth',2); hold on;
xlabel('$\bar{y}$','interpreter','latex')
ylim([0,0.08])
box off
s4.FontSize = fs;
s4.TickDir='out'; 
s4.YTickLabel = '';
s4.Position = [0.74,0.2,0.18,0.8];
fh.PaperPosition=[0 0 16 4];
print('-depsc',['./figures/distributions'])