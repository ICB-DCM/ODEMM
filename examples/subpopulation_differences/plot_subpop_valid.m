% Visualization script for the validation Figure 5D-F of Loos et al., Cell
% Systems (2018).

close all
clear all
clc

load_plot_settings
%% TrkA levels
load('./data/data_PDL');
for r = 1:4
    TP = squeeze(D(3).replicate{r}(end,1,:,:));
    T = TP(:,1);
    P = TP(:,2);
    [Ps,ind] = sort(P);
    Ts = T(ind);
    PTthres = 2600;
    indresp = find(Ps>PTthres,1,'first');
    Tnon(r) = nanmean(Ts(1:indresp));
    Tres(r) = nanmean(Ts(indresp+1:end));
end
b=bar([mean(Tnon),mean(Tres)]);
b.FaceColor = color.param(5,:); hold on;
myerrorbar(1,mean(Tnon),std(Tnon),0.1);
myerrorbar(2,mean(Tres),std(Tres),0.1);
[hy,pv] = ttest2(Tnon,Tres,'Vartype','Unequal')
[hy2,pv2] = ttest(Tnon-Tres)

box off
xlim([0.5,2.5])
set(gca,'xticklabel',{'pErk-','pErk+'},'TickDir','out','Fontsize',fs)
ylabel('TrkA levels [UI]','Fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2.5 2])
%print('-dpdf','./figures/subpop_diff_TrkA')
%% Erk1/2 expression 
close all
for r = 1:4
    EP = squeeze(D(4).replicate{r}(end,1,:,:));
    E = EP(:,1);
    P = EP(:,2);
    [Ps,ind] = sort(P);
    Es = E(ind);
    PEthres = 2600;
    indresp = find(Ps>PEthres,1,'first');
    Enon(r) = nanmean(Es(1:indresp));
    Eres(r) = nanmean(Es(indresp+1:end));
end
b=bar([mean(Enon),mean(Eres)]);
b.FaceColor = color.param(6,:); hold on;
myerrorbar(1,mean(Enon),std(Enon),0.1);
myerrorbar(2,mean(Eres),std(Eres),0.1);
[hy,pv] = ttest2(Enon,Eres,'Vartype','Unequal')
[hy2,pv2] = ttest(Enon-Eres)

box off
xlim([0.5,2.5])
set(gca,'xticklabel',{'pErk-','pErk+'},'TickDir','out','Fontsize',fs)
ylabel('Erk1/2 levels [UI]','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2.5 2])
%print('-dpdf','./figures/subpop_diff_Erk')
%% Erk1/2 dephosphorylation
% for the analysis see ./dephosphorylation/main_analyze_dephospho.m
close all
load ./dephosphorylation/results_dephospho
b=bar([mean(k5_R),mean(k5_T)]);
b.FaceColor = color.param(4,:); hold on;
myerrorbar(1,mean(k5_R),std(k5_R),0.1);
myerrorbar(2,mean(k5_T),std(k5_T),0.1);
[hy,pv] = ttest2(k5_R,k5_T,'Vartype','Unequal')
box off
xlim([0.5,2.5])
set(gca,'xticklabel',{'pErk-','pErk+'},'TickDir','out','Fontsize',fs)
ylabel('Erk1/2 dephosphorylation [min^{-1}]','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2.5 2])
%print('-dpdf','./figures/subpop_diff_k5')




