% Visualization script for ranking of the models and differences

clear all
close all
clc

load_plot_settings

load('./project/data/data_matrices_1D2D');
load ./project/results/results_4_1802

n_data = 0;
for e = 1:8
    n_data = n_data + sum(sum(sum(~isnan(D(e).y))));
end

logPosts = nan(1,128);
BICs = nan(1,128);
AICs = nan(1,128);
conv = nan(1,128);
started = nan(1,128);

notfinished = [];
for i = 1:128
    try
        if isempty(parameters{i})
            notfinished = [notfinished,i];
        else
            logPosts(i) = parameters{i}.MS.logPost(1);
            parameters{i}.MS.AIC = -2*parameters{i}.MS.logPost + 2*parameters{i}.number;
            AICs(i) =  parameters{i}.MS.AIC(1);
            parameters{i}.MS.BIC = -2*parameters{i}.MS.logPost + log(n_data)*parameters{i}.number;
            BICs(i) =  parameters{i}.MS.BIC(1);
            conv(i) = sum(2*(parameters{i}.MS.logPost-parameters{i}.MS.logPost(1))>-icdf('chi2',0.95,1));
            started(i) = sum(~isnan(parameters{i}.MS.logPost));
            tcpus(i) = nanmean(parameters{i}.MS.t_cpu);
        end
    catch
        notfinished = [notfinished,i];
    end
end

comb = nan(2^7,7);
i = 1;
for i1 = 0:1
    for i2 = 0:1
        for i3 = 0:1
            for i4 = 0:1
                for i5 = 0:1
                    for i6 = 0:1
                        for i7 = 0:1
                            comb(i,:) = [i1,i2,i3,i4,i5,i6,i7];
                            i = i+1;
                        end
                    end
                end
            end
        end
    end
end
[BICsort,inds] = sort(BICs);
comb_temp = bsxfun(@times,comb(inds,:),[6:-1:1,7]);
comb_temp(comb_temp == 0) = NaN;

%% BIC weights
figure('name','BIC weights');
sortcomb = comb(inds,:);
diffs_b = BICsort(1:end)-BICsort(1);
ws = exp(-0.5*diffs_b)./nansum(exp(-0.5*diffs_b));
c = sortcomb;
probs = nansum(bsxfun(@times,ws,c'),2);
[probs,ranking] = sort(probs,'descend');

b = bar([probs';ones(1,7)]);
for i = 1:7
    b(i).FaceColor = color.param(ranking(i),1:3);
end
xlim([0.6,1.4])
ylim([0,1]);
box on
set(gca,'xticklabel','')
set(gca,'Position',[0.2,0.2,0.7,0.7]);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 3.9 1.5])
print('-depsc',['./project/figures/ECM_T_BICweights'])


%%
figure('name','ranked models')
rankval = nan(1,7);
for i = 1:7
    temp = BICsort'.*~isnan(comb_temp(:,i));
    rankval(i) = nanmean([1:128]'.*comb_temp(:,i)/nanmean(comb_temp(:,i)));
end
for i = 1:7
    temp = comb_temp(:,ranking(i));
    temp(~isnan(temp)) = 8-i;
    plot(1:128,temp,'.','MarkerSize',8,'Color',color.param(ranking(i),1:3)); hold on;
    h= plot(rankval(ranking(i)),8-i,'ko','MarkerSize',4);
end
xlim([1,128])
ylim([0.5,7.5])
set(gca,'ytick',[1:8]);
set(gca,'yticklabel',{'w','k_4','k_1','k_2','k_3[TrkA]_0','k_5','c[Erk]_0'});
set(gca,'xtick',[1,16,32,48,64,80,96,112,128]);
ylabel({'parameter differences', 'between scaffolds'},'FontSize',6)
xlabel('model ranking according to BIC','fontsize',6)
grid on
set(gca,'xticklabel',{'1','','','','64','','','','128'},'xlim',[1,128]);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7.5 3.5])
print('-dpdf','./project/figures/ECM_rankval')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 14 3])
print('-dpdf','./project/figures/ECM_rankval_Suppl')
%% BIC values
figure('name','BICs');
[BICsort] = sort(BICs);
plot(1:128,BICsort-min(BICsort),'.k','MarkerSize',8); hold on;
set(gca,'xtick',[1,16,32,48,64,80,96,112,128],'xlim',[1,128],'ylim',[-100,2500]);
ylabel('BIC-min(BIC)','FontSize',6)
xlabel('ranked models','fontsize',6)
grid on
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 14 3])
print('-dpdf','./project/figures/BICs')
%% Zoom in for BIC
figure('name','BICs zoom');
plot(1:16,BICsort(1:16)-min(BICsort),'.k','MarkerSize',8); hold on;
set(gca,'xtick',[1,8,16],'xlim',[1,16])%,'ylim',[-10,150]);
ylabel('BIC-min(BIC)','FontSize',6)
xlabel('ranked models','fontsize',6)
grid off
plot([1,16],[10,10],'k:')
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5 3])
print('-dpdf','./project/figures/BICs_zoom')
