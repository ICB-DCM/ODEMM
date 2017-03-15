% main function for the scaling study for number of data points 
% 100,500,1000,5000,1e4,1e5

clear all
close all
clc
global timeval
timeval = 0;

global timesim
timesim = 0;

parameters.name = {'log_{10}(k_{1,1})','log_{10}(k_{1,2})','log_{10}(k_2)',...
    'log_{10}m(log(k_3))','log_{10}(\sigma_{k_3}))',...
    'log_{10}(\sigma_{noise}))',...
    'w'};
parameters.number = length(parameters.name);
parameters.min = [-3*ones(4,1);-3;-3;0];
parameters.max = [ 3*ones(4,1);2;2;1];
M.name = 'CR_SP_k3';
eval(['ODEMM_' M.name]);
options.simulate_musigma = 1;
options.MS.fmincon = optimset('GradObj','on','display','iter','TolFun',1e-10, 'TolX',1e-10, 'MaxIter', 1000,'algorithm','interior-point');
options.MS.n_starts = 10;
options.MS.comp_type = 'sequential'; options.MS.mode = 'visual';

n_data = [100,1000,1e4,1e5];
for iteration = 1:3
    load(['./project/data/conversionprocess_data_' num2str(n_data(1))])
    [conditions,D] = collectConditions(D,M);
    par0 =   getParameterGuesses(parameters,@(xi) logLikelihood(xi,M,D,options,conditions),...
        options.MS.n_starts, parameters.min,parameters.max);
    parameters.guess = par0;
    for n = 1:numel(n_data)
        load(['./project/data/conversionprocess_data_' num2str(n_data(n))])
        [conditions,D] = collectConditions(D,M);
        timeval = 0;
        timesim = 0;
        param{iteration,n} = getMultiStarts(parameters,@(xi) logLikelihood(xi,M,D,options,conditions),options.MS);
        param{iteration,n}.timeval = timeval;
        param{iteration,n}.timesim= timesim;
        param{iteration,n}.n_data = n_data(n);
        save(['./project/results/cr_scalingstudy_simtime_' num2str(iteration) '_' num2str(n)])
        clear D conditions
    end
end
        save(['./project/results/cr_scalingstudy_simtime'])


%%
load_plot_settings
load('./project/results/cr_scalingstudy_simtime')
close all
for n = 1:numel(n_data)
    for iteration = 1:3       
        evaltime(iteration,n) = param{iteration,n}.timeval;
        simtime(iteration,n) = nansum(param{iteration,n}.timesim);
        tcpus(iteration,n) = nansum(param{iteration,n}.MS.t_cpu);
        iters(iteration,n) = nansum(param{iteration,n}.MS.n_iter);
        objevals(iteration,n) = nansum(param{iteration,n}.MS.n_objfun);

        param{iteration,n}.n_data
    end
    n_data(n) = param{iteration,n}.n_data;
end
figure
%errorbar(n_data,mean(tcpus),std(tcpus))
ind = 1:4;%[1,3,5,6];
loglog(n_data(ind),mean(tcpus(:,ind)),'-o','MarkerSize',4); hold on;
loglog(n_data(ind),mean(simtime(:,ind)),'-o','MarkerSize',4)
loglog(n_data(ind),mean(evaltime(:,ind)),'-o','MarkerSize',4)
%loglog(n_data(ind),mean(evaltime(:,ind))+mean(simtime(:,ind)),'-o','MarkerSize',4)

box off
xlabel('number of cells per time point')
ylabel('cpu time [min]')
set(gca,'fontsize',fs,'tickdir','out')
%leg = legend('total time', 'simulation time','density evaluation');
%set(leg,'fontsize',fs,'location','southeast')
set(gca,'Position',[0.3,0.3,0.6,0.6])
%%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 10 6])
feval('print', '-dpdf',['./project/figures/cr_scaling']);
