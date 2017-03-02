% Visualization script for the validation Figure 5D-F

close all
clear all
clc

load_plot_settings
%% auf replicates
load('./project/data/data_Lys_1D2D_sepscaled_withrepl');
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
%[hy,pv]=ttest(Tnon-Tres)
ttest2(Tnon,Tres)
box off
xlim([0.5,2.5])
set(gca,'xticklabel',{'pErk-','pErk+'},'TickDir','out','Fontsize',fs)
ylabel('TrkA levels [UI]','Fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2.5 2])
print('-dpdf','./project/figures/subpop_diff_TrkA')
%%
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
%[hy,pv]=ttest(Enon-Eres)
ttest2(Enon,Eres)
box off
xlim([0.5,2.5])
set(gca,'xticklabel',{'pErk-','pErk+'},'TickDir','out','Fontsize',fs)
ylabel('Erk1/2 levels [UI]','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2.5 2])
print('-dpdf','./project/figures/subpop_diff_Erk')
%%
close all
load ./dephosphorylation/k5s_onlyk5
b=bar([mean(k5_R),mean(k5_T)]);
b.FaceColor = color.param(4,:); hold on;
myerrorbar(1,mean(k5_R),std(k5_R),0.1);
myerrorbar(2,mean(k5_T),std(k5_T),0.1);
ttest2(k5_R,k5_T)
%[hy,pv]=ttest(Enon-Eres)
box off
xlim([0.5,2.5])
set(gca,'xticklabel',{'pErk-','pErk+'},'TickDir','out','Fontsize',fs)
ylabel('dephosphorylation [min^{-1}]','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2.5 2])
print('-dpdf','./project/figures/subpop_diff_k5')



% load('./project/data/data_Lys_1D2D_sepscaled');
% %% TrkA
% TP = squeeze(D(3).y(end,1,:,:));
% T = TP(:,1);
% P = TP(:,2);
% 
% [Ps,ind] = sort(P);
% Ts = T(ind);
% mT = [];
% PTthres = 2550;
% indresp = find(Ps>PTthres,1,'first');
% 
% grid = logspace(2,log10(Ps(end)),15);
% for c = 1:numel(grid)
%     indth = find(Ps<grid(c));
%     mT(c) = nanmean(Ts(indth));
% end
% % %%
% % figure
% % nanmean(Ts(1:indresp))
% % nanmean(Ts(indresp+1:end))
% % ttest2(Ts(1:indresp),Ts(indresp+1:end))
% %%
% % figure
% % loglog(grid,mT,'Marker','.','MarkerSize',8,'Color',color.param(5,:),'LineWidth',1.2);
% % set(gca,'TickDir','out')
% % set(gca,'xtick',[0,1e2,1e3,1e4])%,'xticklabel',{'0','1','2'})
% % xlabel('pErk threshold','FontSize',fs)
% % ylabel('mean TrkA levels','FontSize',fs)
% % xlim([1e2,2e4])
% % box off
% % set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 3 2.5])
% % print('-dpdf','./project/figures/subpop_diff_valid_TrkA')
% %% Erk
% EP = squeeze(D(4).y(end,1,:,:));
% E = EP(:,1);
% P = EP(:,2);
% [Ps,ind] = sort(P);
% Es = E(ind);
% indresp = find(Ps>2723,1,'first');
% grid = logspace(2,log10(Ps(end)),15);
% mE = [];
% for c = 1:numel(grid)
%     indth = find(Ps<grid(c));
%     mE(c) = nanmean(Es(indth));
% end
% %%
% % figure
% % plot(log(E),log(P),'o');
% %%
% bar([nanmean(Es(1:indresp)),nanmean(Es(indresp+1:end))])
% ttest2(Es(1:indresp),Es(indresp+1:end))
% 
% %%
% figure
% loglog(grid,mE,'Marker','.','MarkerSize',8,'Color',color.param(6,:),'LineWidth',1.2);
% set(gca,'TickDir','out')
% set(gca,'xtick',[0,1e2,1e3,1e4]);%,'xticklabel',{'0','1','2'})
% xlabel('pErk threshold','FontSize',fs)
% ylabel('mean Erk levels','FontSize',fs)
% ylim([1e2,1.1*10^3])
% box off
% xlim([1e2,2e4])
% set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 3 2.5])
% print('-dpdf','./project/figures/subpop_diff_valid_Erk')
% 
% %%
% % t = D.t
% % for r = 1:4
% %    % subplot(2,2,r);
% %     for i = 1:numel(t)
% %         y = D.WFinh(r).TrkA{i};
% %         y(y<=0) = NaN;
% %         y2 = D.ctrlinh(r).TrkA{i};
% %         y2(y2<=0) = NaN;
% %         pErk_T(r,i) = nanmean(D.WFinh(r).pErk{i}((y>670)))-nanmean(D.ctrlinh(r).pErk{i}((y2>670)));
% %         pErk_R(r,i) = nanmean(D.WFinh(r).pErk{i}((y<630)))-nanmean(D.ctrlinh(r).pErk{i}((y2<630)));
% %     end
% %
% % %         pErk_T(r,:) = pErk_T(r,:)-min(pErk_T(r,:));
% % %         pErk_T(r,:) =  pErk_T(r,:)./pErk_T(r,1);
% % %         pErk_R(r,:) = pErk_R(r,:)-min(pErk_R(r,:));
% % %         pErk_R(r,:) =  pErk_R(r,:)./pErk_R(r,1);
% %         pErk_T(r,:) = pErk_T(r,:)-(pErk_T(r,end));
% %         pErk_T(r,:) =  pErk_T(r,:)./pErk_T(r,1);
% %         pErk_R(r,:) = pErk_R(r,:)-(pErk_R(r,end));
% %         pErk_R(r,:) =  pErk_R(r,:)./pErk_R(r,1);
% % %      plot(t,pErk_T(r,:),'r'); hold on;
% % %     plot(t,pErk_R(r,:),'b'); hold on;
% % %     xlim([0,37]);
% % %     ylim([-0.2,1])
% % %     set(gca,'ytick',[0,1])
% % %     title(['\rm replicate ' num2str(r)],'FontSize',6);
% % %     if r < 3
% % %     set(gca,'xticklabel','')
% % %     end
% % %     box off
% % %     set(gca,'TickDir','out')
% % end
% % figure
% % errorbar(t, nanmean(pErk_T), std(pErk_T),'Color',[87,104,51]./255,'LineWidth',0.8); hold on;
% % %plot(t, nanmean(pErk_T), '.','MarkerSize',5,'Color',[87,104,51]./255); hold on;
% % %plot(t, nanmean(pErk_T),'Color',[87,104,51]./255); hold on;
% %
% % errorbar(t+0.2, nanmean(pErk_R), std(pErk_R),'Color',color.param(4,:),'LineWidth',0.8); hold on;
% % %plot(t, nanmean(pErk_R), '.','MarkerSize',5,'Color',color.param(4,:)); hold on;
% % %plot(t, nanmean(pErk_R),'Color',color.param(4,:)); hold on;
% %
% % set(gca,'TickDir','out')
% % set(gca,'xtick',[0,20,40]);%,'xticklabel',{'0','1','2'})
% % xlabel('time [min]','FontSize',fs)
% % ylabel({'scaled mean', 'Erk levels [UI]'},'FontSize',fs)
% % %ylim([1e2,1.1*10^3])
% % box off
% % xlim([0,40])
% % ylim([-0.2,1.05])
% % %legend('TrkA+','Ret+');
% % set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2 2])
% % print('-dpdf','./project/figures/subpop_diff_valid_k5')


