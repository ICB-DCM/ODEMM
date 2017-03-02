% Visualization script for the fits in the supplement

close all
clear all
clc

load('./project/data/data_matrices_1D2D');

load ./project/results/results_4_1802
load_plot_settings
p = parameters{15};
clear parameters;
parameters = p;
ODEMM_NGF_ECM_T_comb15
%%

TextSizes.DefaultAxesFontSize = 6;
TextSizes.DefaultTextFontSize = 6;
set(0,TextSizes);

options_plot.x_scale = 'log';
options_plot.data.plot = 'filled';
options_plot.data.bins = 40;
options_plot.model.points = 200;
options_plot.data.lw = 1;
options_plot.model.lw = 0.8;


options_plot.data.col{1} = color.Lys_data;
options_plot.model.col{1} =  color.Lys;
options_plot.model.col{3} =  color.Lys;
options_plot.model.col{5} =  color.Lys;
options_plot.model.col{7} =  color.Lys;

options_plot.data.col{2} = color.ColI_data;
options_plot.model.col{2} = color.ColI;
options_plot.model.col{4} =  color.ColI;
options_plot.model.col{6} =  color.ColI;
options_plot.model.col{8} =  color.ColI;

options_plot.data.col{1} = color.Lys_data;
options_plot.model.col{1} =  color.Lys;
options_plot.model.col{3} =  color.Lys;
options_plot.model.col{5} =  color.Lys;
options_plot.model.col{7} =  color.Lys;

options_plot.data.col{2} = color.ColI_data;
options_plot.model.col{2} = color.ColI;
options_plot.model.col{4} =  color.ColI;
options_plot.model.col{6} =  color.ColI;
options_plot.model.col{8} =  color.ColI;
options_plot.data.col{3} = options_plot.data.col{1};
options_plot.data.col{5} = options_plot.data.col{1};
options_plot.data.col{7} = options_plot.data.col{1};
options_plot.data.col{4} = options_plot.data.col{2};
options_plot.data.col{6} = options_plot.data.col{2};
options_plot.data.col{8} = options_plot.data.col{2};

options_plot.data.fill_col{1} = options_plot.data.col{1};
options_plot.data.fill_col{2} = options_plot.data.col{2};
options_plot.data.fill_col{3} = options_plot.data.col{3};
options_plot.data.fill_col{4} = options_plot.data.col{4};

for d = 1:6
    options_plot.model.levelsets{5,d} = [.45 0.75 1.15 1.5];
    options_plot.model.levelsets{6,d} = [.45 0.75 1.15 1.5];
    
    options_plot.model.levelsets{7,d} = [0.9,1.7,2.5,3.3];
    options_plot.model.levelsets{8,d} = [0.9,1.7,2.5,3.3];
    
end
options_plot.model.levelsets{7,1} = [1.3,2.6,4,5];
options_plot.model.levelsets{8,1} = [1.3,2.6,4,5];
options_plot.fit_scaling_parameters = true;
options_plot.boundaries(1).y_min = 1e2; options_plot.boundaries(1).y_max = 2e4;
options_plot.boundaries(2).y_min = 1e2; options_plot.boundaries(2).y_max = 2e4;
options_plot.boundaries(3).y_min = 1e2; options_plot.boundaries(3).y_max = 2e4;
options_plot.boundaries(4).y_min = 1e2; options_plot.boundaries(4).y_max = 2e4;

options_plot.boundaries(5).y_min = [1e1,5e1];
options_plot.boundaries(5).y_max = [2e4,2e4];
options_plot.boundaries(6).y_min = [1e1,5e1];
options_plot.boundaries(6).y_max = [2e4,2e4];
options_plot.boundaries(7).y_min = [1e2,1e2];
options_plot.boundaries(7).y_max = [5e3,5e4];
options_plot.boundaries(8).y_min = [1e2,1e2];
options_plot.boundaries(8).y_max = [5e3,5e4];
options_plot.subplot_lin = 1;

options_plot.switch_axes = true;
options_plot.marginals = 0;
options_plot.sameplot = 0;
options_plot.model.subpopulations = 0;
options_plot.simulate_musigma = true;
options_plot.data.markersize = 10;
options_plot.model.level_linewidth = 0.5;
xi = parameters.MS.par(:,1);
[conditions,D] = collectConditions(D,M);
close all
[fh] = plotODEMix(D,M,xi,[5:8],options_plot);
%[fh] = plotODEMix(D,[],[],[3:4],options_plot,tu_ind);

mtvals(1:8) = log10([2:9]*10^1);
mtvals(9:16) = log10([2:9]*10^2);
mtvals(17:24) = log10([2:9]*10^3);
mtvals(25:32) = log10([2:9]*10^4);

%% TrkA/pErk
figure(fh{5});

for i = 12:-1:1
    s=subplot(2,6,i)
    view(180,90);
    set(gca,'xdir','reverse');
    if i <= 6
        set(gca, 'xtick', [2,3,4] ,'xticklabel',{''},'ytick', [2,3,4],...
            'yticklabel',{''});
        
        hA = gca;
        hA.YAxis.MinorTickValues = mtvals;
        hA.XAxis.MinorTickValues = mtvals;
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(s, 'position', [0.1+(i-1)*0.14,0.5,0.13,0.25] );
        
    elseif i > 6
        set(gca, 'xtick', [1e2,1e3,1e4] ,'ytick',[1e2,1e3,1e4],'yticklabel',{''},...
            'xticklabel',{'10^2','','10^4'});
        set(s, 'position', [0.1+(i-7)*0.14,0.2,0.13,0.25] );
        
    end
    
    if i == 1
        set(gca  ,'yticklabel',{'10^2','10^3','10^4'});
    elseif i == 7
        set(gca,'yticklabel',{'10^2','10^3','10^4'});
    end
end
%%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 11 5.5])
set(gcf,'renderer','opengl');
print('-dpdf','-r500','./project/figures/datafit_TrkApErk_PDL')

figure(fh{6});

for i = 12:-1:1
    s=subplot(2,6,i)
    view(180,90);
    set(gca,'xdir','reverse');
    if i <= 6
        set(gca, 'xtick', [2,3,4] ,'xticklabel',{''},'ytick', [2,3,4],...
            'yticklabel',{''});
        
        hA = gca;
        hA.YAxis.MinorTickValues = mtvals;
        hA.XAxis.MinorTickValues = mtvals;
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(s, 'position', [0.1+(i-1)*0.14,0.5,0.13,0.25] );
        
    elseif i > 6
        set(gca, 'xtick', [1e2,1e3,1e4] ,'ytick',[1e2,1e3,1e4],'yticklabel',{''},...
            'xticklabel',{'10^2','','10^4'});
        set(s, 'position', [0.1+(i-7)*0.14,0.2,0.13,0.25] );
        
    end
    
    if i == 1
        set(gca  ,'yticklabel',{'10^2','10^3','10^4'});
    elseif i == 7
        set(gca,'yticklabel',{'10^2','10^3','10^4'});
    end
end
%%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 11 5.5])
set(gcf,'renderer','opengl');
print('-dpdf','-r500','./project/figures/datafit_TrkApErk_ColI')
%% Erk/pErk
figure(fh{7});
for i = 12:-1:1
    s=subplot(2,6,i)
    view(180,90);
    set(gca,'xdir','reverse');
    if i <= 6
        set(gca, 'xtick', [2,3] ,'xticklabel',{''},'ytick', [2,3,4],...
            'yticklabel',{''});
        hA = gca;
        hA.YAxis.MinorTickValues = mtvals;
        hA.XAxis.MinorTickValues = mtvals;
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(s, 'position', [0.1+(i-1)*0.14,0.5,0.13,0.25] );
        
    elseif i > 6
        set(gca, 'xtick', [1e2,1e3] ,'ytick',[1e2,1e3,1e4],'yticklabel',{''},...
            'xticklabel',{'10^2','10^3'});
        set(s, 'position', [0.1+(i-7)*0.14,0.2,0.13,0.25] );
        
    end
    
    if i == 1
        set(gca  ,'yticklabel',{'10^2','10^3','10^4'});
    elseif i == 7
        set(gca,'yticklabel',{'10^2','10^3','10^4'});
    end
end
%%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 11 5.5])
set(gcf,'renderer','opengl');
print('-dpdf','-r500','./project/figures/datafit_ErkpErk_PDL')

%% Erk/pErk
figure(fh{8});
for i = 12:-1:1
    s=subplot(2,6,i)
    view(180,90);
    set(gca,'xdir','reverse');
    if i <= 6
        set(gca, 'xtick', [2,3] ,'xticklabel',{''},'ytick', [2,3,4],...
            'yticklabel',{''});
        hA = gca;
        hA.YAxis.MinorTickValues = mtvals;
        hA.XAxis.MinorTickValues = mtvals;
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(s, 'position', [0.1+(i-1)*0.14,0.5,0.13,0.25] );
        
    elseif i > 6
        set(gca, 'xtick', [1e2,1e3] ,'ytick',[1e2,1e3,1e4],'yticklabel',{''},...
            'xticklabel',{'10^2','10^3'});
        set(s, 'position', [0.1+(i-7)*0.14,0.2,0.13,0.25] );
        
    end
    
    if i == 1
        set(gca  ,'yticklabel',{'10^2','10^3','10^4'});
    elseif i == 7
        set(gca,'yticklabel',{'10^2','10^3','10^4'});
    end
end
%%
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 11 5.5])
set(gcf,'renderer','opengl');
print('-dpdf','-r500','./project/figures/datafit_ErkpErk_ColI')
close all

%%
close all
tu_ind{1} = [1:7];
tu_ind{2} = [1:7];
tu_ind{3} = [1:7];
tu_ind{4} = [1:7];
options_plot.sameplot = 1;
fh = plotODEMix_ECM(D,M,xi,[1,2],options_plot,tu_ind);
figure(fh{1})
for i = 7:-1:1
    s = subplot(1,7,i);
    set(gca,'ylim',[0,0.15])   
end
set(gca,'TickDir','out')
ah=axes('position',[.13,.11,0.775,0],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.5,'Color','k');
ah=axes('position',[.13,.11,0,0.8],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.5,'Color','k');
figure(fh{1})
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.7 1.7])
set(gcf,'renderer','opengl');
feval('print', '-dpdf','-r1000','./project/figures/datafitPkin_suppl.pdf');
%%
close all
fh = plotODEMix_ECM(D,M,xi,[3,4],options_plot,tu_ind);
figure(fh{3})
for i = 7:-1:1
    s = subplot(1,7,i);
    set(gca,'ylim',[0,0.15])   
end
set(gca,'TickDir','out')
ah=axes('position',[.13,.11,0.775,0],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.5,'Color','k');
ah=axes('position',[.13,.11,0,0.8],'visible','off');
line([0,1],[0,1],'parent',ah,'linewidth',0.5,'Color','k');
figure(fh{3})
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.7 1.7])
set(gcf,'renderer','opengl');
feval('print', '-dpdf','-r1000','./project/figures/datafitPdr_suppl.pdf');

