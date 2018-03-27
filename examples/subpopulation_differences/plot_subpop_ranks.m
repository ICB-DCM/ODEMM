% Visualization script for the ranking of the 64 models accounting for
% differences between subpopulations.

clear all
close all
clc

load ./data/data_PDL
load ./results/results_subpop_differences
load_plot_settings
n_data = 0;
for e = 1:4
    n_data = n_data + sum(sum(sum(~isnan(D(e).y))));
end
AICs = nan(1,64);
logPosts = nan(1,64);
BICs = nan(1,64);
converged = nan(1,64);
difference_names = {'k_1','k_2','k_4','k_5','k_3[TrkA]_0','c[Erk]_0'};
ranking = [5,6,4,3,1,2];
for i = 1:64
    try
        logPosts(i) = parameters{i}.MS.logPost(1);
        parameters{i}.MS.AIC = -2*parameters{i}.MS.logPost + 2*parameters{i}.number;
        AICs(i) =  parameters{i}.MS.AIC(1);
        parameters{i}.MS.BIC = -2*parameters{i}.MS.logPost + log(n_data)*parameters{i}.number;
        BICs(i) =  parameters{i}.MS.BIC(1);
        converged(i) = sum(2*(parameters{i}.MS.logPost-parameters{i}.MS.logPost(1))>-...
            icdf('chi2',0.95,1));
    catch
    end
end
find(converged>0 & converged<3)

%% combinations
comb = nan(2^6,6);
i = 1;
for i1 = 0:1
    for i2 = 0:1
        for i3 = 0:1
            for i4 = 0:1
                for i5 = 0:1
                    for i6 = 0:1
                        comb(i,:) = [i1,i2,i3,i4,i5,i6];
                        i = i+1;
                    end
                end
            end
        end
    end
end
%% sort BIC values
figure('name','ranked importance');
[BICsort,inds] = sort(BICs);
comb_temp = bsxfun(@times,comb(inds,:),6:-1:1);
comb_temp(comb_temp == 0) = NaN;
rankval = nan(6,1);
for i = 1:6
    rankval(i) = nanmean([1:64]'.*comb_temp(:,i)/nanmean(comb_temp(:,i)));
end
for i = 1:6
    temp = comb_temp(:,ranking(i));
    temp(~isnan(temp)) = 7-i;
    plot(1:64,temp,'.','MarkerSize',8,'Color',color.param(ranking(i),1:3)); hold on;
    h= plot(rankval(ranking(i)),7-i,'ko','MarkerSize',3,'LineWidth',1);
end
for i = 1:6
    yticklabelnames{7-i} = difference_names{ranking(i)};
end
set(gca,'xtick',[1,16,32,48,64],'ytick',[1:6],'yticklabels',yticklabelnames,...
    'xlim',[1,64],...
    'ylim',[0.5,6.5] );
ylabel({'subpopulation difference'},'FontSize',6)
xlabel('model ranking according to BIC','fontsize',6)
grid on
legend(h,'mean rank')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7.5 2.5])
%print('-dpdf','./figures/subpop_rankval')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 14.05 2.5])
%print('-dpdf','./figures/subpop_rankval_Suppl')
%% BIC values
figure('name','BICs');
[BICsort,inds] = sort(BICs);
plot(1:64,BICsort-min(BICsort),'.k','MarkerSize',8); hold on;
set(gca,'xtick',[1,16,32,48,64],'xlim',[1,64],'ylim',[-500,2e4]);
ylabel('BIC-min(BIC)','FontSize',6)
xlabel('ranked models','fontsize',6)
grid on
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2 2.5])
%print('-dpdf','./figures/BICs')
%% Zoom for BIC values
figure('name','BICs_zoom');
[BICsort,inds] = sort(BICs);
plot(1:16,BICsort(1:16)-min(BICsort),'.k','MarkerSize',8); hold on;
set(gca,'xtick',[1,8,16],'xlim',[1,16],'ylim',[-10,150]);
ylabel('BIC-min(BIC)','FontSize',6)
xlabel('ranked models','fontsize',6)
grid off
plot([1,16],[10,10],'k:')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5 3])
%print('-dpdf','./figures/BICs_zoom')

%% BIC weights
figure('name','BIC weights');
sortcomb = comb(inds,:);
diffs_b = BICsort-BICsort(1);
ws = exp(-0.5*diffs_b);
ws = ws/nansum(ws);
c = sortcomb;
probs = nansum(bsxfun(@times,ws,c'),2);
probs = probs(ranking);
b = bar([probs';ones(1,6)]);
for i = 1:6
    b(i).FaceColor = color.param(ranking(i),1:3);
end
box on
set(gca,'Position',[0.2,0.2,0.7,0.7],...
    'ytick',[0,0.5,1],...
    'xtick',[],...
    'xticklabel','',...
    'xlim',[0.6,1.4],...
    'ylim',[0,1.1]...
    );

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2.6 1.5])
%print('-depsc','./figures/subpop_BICweights')

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 2.5 1.2])
%print('-depsc','./figures/subpop_BICweights_Suppl')

%% Figure Sampling
figure('name','Sampling');
subplot('Position',[0.1,0.2,0.7,0.1]);
ind_models = [1,2,3,5,9,17,33];
BIConediff = BICs(ind_models);

[sortBIConediff,indsort]=sort(BIConediff,'descend');

load ./results/results_lppds_subpop
load ./results/results_logmargs_subpop

lppds = lppds(indsort);
logmargs = logmargs(indsort);
for m = 1:numel(ind_models)
    plot(m,-lppds(m)+max(lppds),'d',...
    'MarkerSize',4.5,'LineWidth',1.1,'Color',color.samp); hold on;
    plot(m,-logmargs(m)+max(logmargs),'.',...
    'MarkerSize',17,'Color',color.samp); hold on;
    plot(m,0.5*(sortBIConediff(m)-min(sortBIConediff)),'.',...
    'MarkerSize',10,'Color',color.opt); hold on;
    
end
differences = {'none','c[Erk]_0','k_3[TrkA]_0','k_5','k_4','k_2','k_1'};
set(gca,'xticklabel',differences(indsort))
set(gca,'TickDir','out','xlim',[0.9,7.1],'xtick',[1:7],...
    'xlim',[0.5,7.5],'ylim',[-0.1,0.1],'ytick',0);
box off
subplot('Position',[0.1,0.32,0.7,0.5]);
for m = 1:numel(ind_models)
    plot(m,-lppds(m)+max(lppds),'d',...
    'MarkerSize',4.5,'LineWidth',1.1,'Color',color.samp); hold on;
    plot(m,-logmargs(m)+max(logmargs),'.',...
    'MarkerSize',17,'Color',color.samp); hold on;
    plot(m,0.5*(sortBIConediff(m)-min(sortBIConediff)),'.',...
    'MarkerSize',10,'Color',color.opt); hold on;
end

disp(['correlation lppd, BIC: ' num2str(spear(-lppds'+max(lppds),...
   sortBIConediff'-min(sortBIConediff)))])
disp(['correlation BF, BIC: ' num2str(spear(-logmargs'+max(logmargs),...
   sortBIConediff'-min(sortBIConediff)))])
disp(['correlation lppd, BIC: ' num2str(corr(-lppds'+max(lppds), ...
    sortBIConediff'-min(sortBIConediff), 'type', 'Spearman'))])

set(gca,'TickDir','out','yscale','log','xtick',[1:7],...
    'xlim',[0.5,7.5],'yminortick','on','ylim',[1e3,1e4]);
box off
set(gca,'xtick',[],'xticklabel','')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 8 5])
%print('-depsc','./figures/subpop_Sampling_SI')