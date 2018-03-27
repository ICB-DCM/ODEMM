% Main script for the analysis of the dephosphorylation data.

close all
clear all
clc

load data_dephospho
t = D.t;

% In this script, the ending _T stands for the responsive subpoulation and
% the ending _R for the non-responsive subpopulation.

figure(1)
for r = 1:4
    subplot(2,2,r);
    for i = 1:numel(t)
        y = D.WFinh(r).TrkA{i};
        y(y<=0) = NaN;
        y2 = D.ctrlinh(r).TrkA{i};
        y2(y2<=0) = NaN;
        pErk_T(r,i) = nanmean(D.WFinh(r).pErk{i}((y>670)))-nanmean(D.ctrlinh(r).pErk{i}((y2>670)));
        pErk_R(r,i) = nanmean(D.WFinh(r).pErk{i}((y<630)))-nanmean(D.ctrlinh(r).pErk{i}((y2<630)));
    end
    pErk_T(r,:) = pErk_T(r,:)-(pErk_T(r,end));
    pErk_T(r,:) =  pErk_T(r,:)./pErk_T(r,1);
    pErk_R(r,:) = pErk_R(r,:)-(pErk_R(r,end));
    pErk_R(r,:) =  pErk_R(r,:)./pErk_R(r,1);
    plot(t,pErk_T(r,:),'o-r'); hold on;
    plot(t,pErk_R(r,:),'diamond-b'); hold on;
    xlim([0,37]);
    ylim([-0.2,1.3])
    set(gca,'ytick',[0,1])
    title(['\rm replicate ' num2str(r)],'FontSize',6);
    if r < 3
        set(gca,'xticklabel','')
    end
    box off
    set(gca,'TickDir','out')
end

%% Estimate parameters
parameters.name = {'k5'};
parameters.number = 1;
parameters.max = 1;
parameters.min = 0;

options.MS.fmincon = optimset('GradObj','off','display','iter','TolFun',1e-10, 'TolX',1e-14, 'MaxIter', 1000,'algorithm','interior-point');
options.MS.mode = 'text';
options.MS.n_starts = 10;
t = [0,1,4,7,10,13,16,19,22,25,28,31,34,37];

for r = 1:4
    param_R{r} = getMultiStarts(parameters,@(theta) -exp_decay(theta,t,pErk_R(r,:)),options.MS);
    k5_R(r) = param_R{r}.MS.par(1);
    param_T{r} = getMultiStarts(parameters,@(theta) -exp_decay(theta,t,pErk_T(r,:)),options.MS);
    k5_T(r) = param_T{r}.MS.par(1);
end

save results_dephospho k5_R k5_T param_R param_T


