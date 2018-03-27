% Visualization script for the validation of the TrkA and Erk differences

close all
clear all
clc

load_plot_settings

load('./data/data_PDL_ColI'); 
%% TrkA levels
for d = 1:numel(D(5).u)
    for r = 1:4
        TP1 = squeeze(D(5).replicate{r}(d,1,:,:));
        T1 = TP1(:,1);
        
        TP2 = squeeze(D(6).replicate{r}(d,1,:,:));
        T2 = TP2(:,1);
                
        Ts1(r,d) = nanmean(T1);
        Ts2(r,d) = nanmean(T2);
    end
    %[signif(d),pv(d)] = ttest2(Ts1(:,d),Ts2(:,d));
end
[signif_all,pv_all] = ttest2(Ts1(:),Ts2(:),'VarType','Unequal');
[signif_all2,pv_all2] = ttest(Ts1(:)-Ts2(:));

% plot TrkA differences
figure(1)
b=bar([mean(Ts1(:)),mean(Ts2(:));0,0]);
b(1).FaceColor = color.PDL_data; hold on;
b(2).FaceColor = color.ColI_data; hold on;

myerrorbar(0.86,mean(Ts1(:)),std(Ts1(:)),0.05);
myerrorbar(1.145,mean(Ts2(:)),std(Ts2(:)),0.05);
xlim([0.7,1.3])
ylabel('TrkA levels [UI]','fontsize',fs)
set(gca,'xtick',[],'tickdir','out','fontsize',fs)

box off
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2.5 1.8])
%print('-dpdf','./figures/ECM_diff_TrkA')

%% Erk1/2 levels
for d = 1:numel(D(5).u)    
    for r = 1:4
        EP1 = squeeze(D(7).replicate{r}(d,1,:,:));
        E1 = EP1(:,1);
        
        EP2 = squeeze(D(8).replicate{r}(d,1,:,:));
        E2 = EP2(:,1);
        
        Es1(r,d) = nanmean(E1);
        Es2(r,d) = nanmean(E2);
    end
end
[signif_all,pv_all] = ttest2(Es1(:),Es2(:),'Vartype','Unequal');
[signif_all2,pv_all2] = ttest(Es1(:)-Es2(:));

% plot Erk differences
figure(2)
b=bar([mean(Es1(:)),mean(Es2(:));0,0]);
b(1).FaceColor = color.PDL_data; hold on;
b(2).FaceColor = color.ColI_data; hold on;
myerrorbar(0.86,mean(Es1(:)),std(Es1(:)),0.05);
myerrorbar(1.145,mean(Es2(:)),std(Es2(:)),0.05);
xlim([0.7,1.3])
ylabel('Erk1/2 levels [UI]','fontsize',fs)
set(gca,'xtick',[],'tickdir','out','fontsize',fs)
box off
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2.5 1.8])
%print('-dpdf','./figures/ECM_diff_Erk')
