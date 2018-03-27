clear all
close all
clc

load ./data/data_PDL
load('./results/results_subpop_TrkA','parameters','conditions');
load_plot_settings
ODEMM_NGF_subpop_TrkA;

%% initialize estruct for SP
estruct.beta = @(xi) [xi(1);
                  xi(2);
                  xi(3);
                  xi(4);
                  xi(5);
                  xi(6)];
estruct.dbetadxi = @(xi) [1,0,0,0,0,0,0,0,0;
            0,1,0,0,0,0,0,0,0;
            0,0,1,0,0,0,0,0,0;
            0,0,0,1,0,0,0,0,0;
            0,0,0,0,1,0,0,0,0;
            0,0,0,0,0,1,0,0,0];
estruct.delta = @(xi) [xi(7),xi(8),xi(9)];
estruct.ddeltadxi = @(xi)[0,0,0,0,0,0,1,0,0;
             0,0,0,0,0,0,0,1,0;...
             0,0,0,0,0,0,0,0,1];
A = [eye(6)];
B = [zeros(4,2);
    eye(2)];
estruct.phi = @(beta,b) A*beta + B*b;
estruct.dphidbeta = @(beta,b) A;
estruct.dphidb = @(beta,b) B;
estruct.sigma_noise = @(phi) zeros(1,3);
estruct.dsigma_noisedphi = @(phi) zeros(1,3,6);
op_SP.nderiv = 0;
op_SP.req = [1,1,0,0,0];
op_SP.type_D = 'matrix-logarithm';

%% Sampling
w = parameters.MS.par(29,1);
xi = parameters.MS.par(:,1);
tsim = [0,120];
numCells = 1e3;
sc_traj = NaN(numel(tsim),3,numCells);
inds1 = [];
inds2 = [];

for iCell = 1:numCells
    iCell
   % sample w
   if rand(1,1) < w
         theta = M.theta(xi,[1,0])';
         inds1 = [inds1,iCell];
   else
         theta = M.theta(xi,[1,1])';
         inds2 = [inds2,iCell];
   end
   % sample theta for Cell
   beta = estruct.beta(theta);
   delta = estruct.delta(theta);
   [D,~,~,~,~,~] = xi2D(delta,op_SP.type_D);
   n_b = size(D,1);
   bsample = mvnrnd(zeros(n_b,1),D);   
   [~,~,~,sc_traj(:,:,iCell)] = simulate_ODEmodel_sPsET_loglog(tsim,estruct.phi(beta,bsample'),1);    
end

%% Visualization
close all
TrkAlevels = squeeze(sc_traj(1,3,:));
meanTrkA = mean(TrkAlevels);
[~,indsort] = sort(TrkAlevels);
col = linspace(min(TrkAlevels),max(TrkAlevels),1000);
scatter(squeeze(log10(exp(sc_traj(1,1,:)))),squeeze(log10(exp(sc_traj(2,1,:)))),10,TrkAlevels,'filled'); hold on;
load mycmap
colormap(mycmap)
axis square
cb = colorbar;
set(cb,'Ticks',[]);
xlabel('pErk1/2 levels at 0min')
ylabel('pErk1/2 levels at 120min')
box off
set(gca,'TickDir','out')
set(gca,'ytick',[2,3,4],'yticklabel',{'10^2','10^3','10^4'});
set(gca,'xtick',[2,3,4],'xticklabel',{'10^2','10^3','10^4'});

ylim([2,4.5])
xlim([2,3.7])

%%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5 5])
set(gcf,'renderer','opengl');
%print('-dpdf','-r500','./project/figures/sc_prediction')




