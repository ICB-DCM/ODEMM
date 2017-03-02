% Visualization script for the dephosphorylation analysis
% @note todo

close all
for r = 1:4
        %figure(1)
    subplot(1,4,r)
    b = -exp(-k5_R(r)*max(t))/(1-exp(-k5_R(r)*max(t)));
    a = 1-b;
    sim_R = a*exp(-k5_R(r) * t) + b;
  
    b = -exp(-k5_T(r)*max(t))/(1-exp(-k5_T(r)*max(t)));
    a = 1-b;
    sim_T = a*exp(-k5_T(r) * t) + b;
    ldR = plot(t,pErk_R(r,:),'diamond','color','b','MarkerSize',3); hold on;
    ldT = plot(t,pErk_T(r,:),'ro','MarkerSize',3); hold on;
    lsR = plot(t,sim_R,'b','LineWidth',1.2); hold on;  
    lsT = plot(t,sim_T,'r--','LineWidth',1.2); hold on;
    ylim([-0.2,1.2])
        set(gca,'ytick',[0,1],'fontsize',6,'tickdir','out');
    xlabel('time [min]','fontsize',6)
    if r == 1
    ylabel('scaled pErk levels','fontsize',6)
    end
    box off
    text(10,0.95,['k_{5,TrkA-} = ' num2str(k5_R(r))],'color','b','fontsize',6);
    text(10,0.75,['k_{5,TrkA+} = ' num2str(k5_T(r))],'color','r','fontsize',6);
    if r == 4
      %  legend([ldR,ldT,lsR,lsT],'data (TrkA-)','data (TrkA+)','fit (TkrA-)','fit (TrkA+)',...
        %    'orientation','horizontal','fontsize',6)
        
    end
        title(['replicate ' num2str(r)],'fontsize',6)

    
end

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 18 3])
print('-dpdf','./project/figures/dephoshpo_new')

%%
load k5s_new
close all

for r = 1:4
        %figure(1)
    subplot(2,2,r)
    title(['replicate ' num2str(r)])
    theta_R = param_R{r}.MS.par(:,1);
    sim_R = theta_R(1) * exp(-theta_R(2) * t) + theta_R(3);
        lsR = plot(t,sim_R,'b'); hold on;
    ldR = plot(t,pErk_R(r,:),'b*','MarkerSize',4); hold on;
%     figure(2)
%         subplot(2,2,r)
    lsT = plot(t,sim_T,'r'); hold on;

   ldT=  plot(t,pErk_T(r,:),'ro','MarkerSize',4); hold on;
    theta_T = param_T{r}.MS.par(:,1);
    sim_T = theta_T(1) * exp(-theta_T(2) * t) + theta_T(3);
    ylim([-0.2,1.2])
    if r > 2
    xlabel('time [min]')
    else
        set(gca,'xticklabel','')
    end
    if r == 1 || r == 3
    ylabel('scaled pErk')
    else
        set(gca,'yticklabel','')
    end
    set(gca,'ytick',[0,1],'fontsize',6,'tickdir','out');
    box off
    text(20,0.95,num2str(theta_R(2)),'color','b','fontsize',6)
    text(20,0.75,num2str(theta_T(2)),'color','r','fontsize',6)
    if r == 1
        %legend([ldT,ldR,lsT,lsR],'data (TrkA+)','data (Ret+)','fit (TkrA+)','fit (Ret+)')
        
    end
    
end

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 8 6])
print('-dpdf','./project/figures/dephoshpo')
