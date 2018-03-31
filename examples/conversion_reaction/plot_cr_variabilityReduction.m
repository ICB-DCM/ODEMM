% clear all
% close all
% clc
% load_plot_settings
% 
% 
%% SP all
subplot('Position',[0.2,0.64,0.7,0.09])
load('./results/results_SP_all');
try
    load('./results/results_BFchains_SP_all');
    thin = ceil(numel(parameters_l{end}.S.burnin+1:...
            size(parameters_l{end}.S.par,2))/100); 
    means = [];
    for iSample = parameters_l{end}.S.burnin+1:thin:size(parameters_l{end}.S.par,2) 
        s=1;
        t_ind = 1:5;
        xi = parameters_l{end}.S.par(:,iSample,1);
        xi(8) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_total=sigma{s,1}(t_ind);
        xi = parameters_l{end}.S.par(:,iSample,1);
        xi(5) = -Inf;
        xi(8) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_nok1=sigma{s,1}(t_ind);
        xi = parameters_l{end}.S.par(:,iSample,1);
        xi(6) = -Inf;
        xi(8) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_nok2=sigma{s,1}(t_ind);
        xi = parameters_l{end}.S.par(:,iSample,1);
        xi(7) = -Inf;
        xi(8) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_nok3=sigma{s,1}(t_ind);
        means = [means;mean(100*[sigma_nok1./sigma_total,sigma_nok2./sigma_total,sigma_nok3./sigma_total])];
    end
    %%
    b=bar(100-[mean(means)],'BarWidth',0.5); hold on;
    b.FaceColor = [0.8,0.8,0.8]; hold on;
    myerrorbar(1,100-mean(means(:,1)),std(means(:,1)),0.1);
    myerrorbar(2,100-mean(means(:,2)),std(means(:,2)),0.1);
    myerrorbar(3,100-mean(means(:,3)),std(means(:,3)),0.1);
    %ylabel('explained variability [%]','FontSize',fs)
    box off
    set(gca,'xticklabel',{''},'TickDir','out','FontSize',fs,'xlim',[0.65,3.35],...
        'ylim',[0,100]);
catch
    disp('results_BFchains_SP_all not found. Skipped.')
end

%% k2k3
subplot('Position',[0.2,0.49,0.7,0.09])
load('./results/results_SP_k2k3');
try
    load('./results/results_BFchains_SP_k2k3');
    thin = ceil(numel(parameters_l{end}.S.burnin+1:...
            size(parameters_l{end}.S.par,2))/100); 
    means = [];
    for iSample = parameters_l{end}.S.burnin+1:thin:size(parameters_l{end}.S.par,2) 
        s=1;
        t_ind = 1:5;
        xi = parameters_l{end}.S.par(:,iSample,1);
        xi(7) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_total=sigma{s,1}(t_ind);
        xi = parameters_l{end}.S.par(:,iSample,1);
        xi(7) = -Inf;
        xi(5) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_nok2=sigma{s,1}(t_ind);
        xi = parameters_l{end}.S.par(:,iSample,1);
        xi(6) = -Inf;
        xi(7) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_nok3=sigma{s,1}(t_ind);
        means = [means;mean(100*[ones(numel(t_ind),1),sigma_nok2./sigma_total,sigma_nok3./sigma_total])];
    end
    b=bar(100-[mean(means)],'BarWidth',0.5); hold on;
    b.FaceColor = [0.8,0.8,0.8]; hold on;
    myerrorbar(1,100-mean(means(:,1)),std(means(:,1)),0.1);
    myerrorbar(2,100-mean(means(:,2)),std(means(:,2)),0.1);
    myerrorbar(3,100-mean(means(:,3)),std(means(:,3)),0.1);
    %ylabel('explained variability [%]','FontSize',fs)
    ylim([0,100])
    box off
    set(gca,'xticklabel',{''},'TickDir','out','FontSize',fs,'xlim',[0.65,3.35],...
        'ylim',[0,100]);
catch
    disp('results_BFchains_SP_k2k3 not found. Skipped.')
end
%% k1k3 
subplot('Position',[0.2,0.34,0.7,0.09])
load('./results/results_SP_k1k3');
try
    load('./results/results_BFchains_SP_k1k3');
    thin = ceil(numel(parameters_l{end}.S.burnin+1:...
            size(parameters_l{end}.S.par,2))/100); 
    means = [];
    for iSample = parameters_l{end}.S.burnin+1:thin:size(parameters_l{end}.S.par,2) 
        s=1;
        t_ind = 1:5;
        xi = parameters_l{end}.S.par(:,iSample,1);
        xi(7) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_total=sigma{s,1}(t_ind);
        xi = parameters_l{end}.S.par(:,iSample,1);
        xi(7) = -Inf;
        xi(5) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_nok1=sigma{s,1}(t_ind);
        xi = parameters_l{end}.S.par(:,iSample,1);
        xi(6) = -Inf;
        xi(7) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_nok3=sigma{s,1}(t_ind);
        means = [means;mean(100*[sigma_nok1./sigma_total,ones(numel(t_ind),1),sigma_nok3./sigma_total])];
    end
    b=bar(100-[mean(means)],'BarWidth',0.5); hold on;
    b.FaceColor = [0.8,0.8,0.8]; hold on;
    myerrorbar(1,100-mean(means(:,1)),std(means(:,1)),0.1);
    myerrorbar(2,100-mean(means(:,2)),std(means(:,2)),0.1);
    myerrorbar(3,100-mean(means(:,3)),std(means(:,3)),0.1);
    ylim([0,100])
    box off
    set(gca,'xticklabel',{''},'TickDir','out','FontSize',fs,'xlim',[0.65,3.35],...
        'ylim',[0,100]);
catch
    disp('results_BFchains_SP_k1k3 not found. Skipped.')
end
%%
subplot('Position',[0.2,0.19,0.7,0.09])
load('./results/results_SP_k3');
try
    load('./results/results_BFchains_SP_k3');
    thin = ceil(numel(parameters_l{end}.S.burnin+1:...
            size(parameters_l{end}.S.par,2))/100); 
    means = [];
    for iSample = parameters_l{end}.S.burnin+1:thin:size(parameters_l{end}.S.par,2) 
        s=1;
        t_ind = 1:5;
        xi = parameters_l{end}.S.par(:,iSample,1);
        xi(6) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_total=sigma{s,1}(t_ind);
        xi(5) = -Inf;
        [mu,sigma] = evaluateModel(xi,M,D,options);
        sigma_nok3=sigma{s,1}(t_ind);
        means = [means;mean(100*[ones(numel(t_ind),1),ones(numel(t_ind),1),sigma_nok3./sigma_total])];
    end
    b=bar(100-[mean(means)],'BarWidth',0.5); hold on;
    b.FaceColor = [0.8,0.8,0.8]; hold on;
    myerrorbar(1,100-mean(means(:,1)),std(means(:,1)),0.1);
    myerrorbar(2,100-mean(means(:,2)),std(means(:,2)),0.1);
    myerrorbar(3,100-mean(means(:,3)),std(means(:,3)),0.1);
    ylim([0,100])
    box off
    set(gca,'xticklabel',{''},'TickDir','out','FontSize',fs,'xlim',[0.65,3.35],...
        'ylim',[0,100]);
catch
    disp('results_BFchains_SP_k3 not found. Skipped.')
end

