function fh = robustSimStudy_boxplots(BICs,MSEs,all_t_cpus,converged,t_cpus)
fh = figure;
labelFlag = true;
load_plot_settings_robust

% BIC
s1 = subplot('Position',[0.4,0.8,0.5,0.15]);
boxplot(bsxfun(@minus,BICs,min(BICs')')+1,'OutlierSize',4,'Symbol','k.','Width',0.8);hold on;
box off
s1.YScale = 'log';
s1.TickDir = 'out';
s1.XTick = [];
s1.YTick = [1,100,1e4,1e6,1e8];
s1.FontSize = 8;
ylim([8.0001e-1,1e5])
s1.YTickLabel = '';
if labelFlag
    ylabel('\Delta BIC + 1')
else
    s1.YTickLabel = '';
end
color_boxplots_robustSimStudy

% Mean Squared Error
s2 = subplot('Position',[0.4,0.62,0.5,0.15]);
boxplot(MSEs,'OutlierSize',4,'Symbol','k.','Width',0.8);hold on;
s2.YScale = 'log';
s2.TickDir = 'out';
s2.XTick = [];
s2.FontSize = 8;
s2.YTick = [1e-4,1e-2,1e0,1e2];
if labelFlag
    ylabel('MSE')
else
    s2.YTickLabel = '';
end

ylim([5e-6,1e2])
color_boxplots_robustSimStudy

% tcpus per run
s3 = subplot('Position',[0.4,0.44,0.5,0.15]);
boxplot([all_t_cpus{1},all_t_cpus{2},all_t_cpus{3},all_t_cpus{4}],...
    'OutlierSize',4,'Symbol','k.','Width',0.8);hold on;
s3.TickDir = 'out';
s3.FontSize = 8;
s3.XTick = [];
s3.YTick = [1e-1,1e0,1e1,1e2];
s3.YScale = 'log';
ylim([6e-2,5e2])
if labelFlag
    ylabel('CPU time per start')
else
    s3.YTickLabel = '';
end
color_boxplots_robustSimStudy

% converged
s4 = subplot('Position',[0.4,0.26,0.5,0.15]);
boxplot(100*converged/30,'OutlierSize',4,'Symbol','k.','Width',0.8);hold on;
ylim([0,100]);
s4.TickDir = 'out';
s4.FontSize = 8;
s4.XTick = [];
if labelFlag
    ylabel('converged starts (%)')
else
    s4.YTickLabel = '';
end
color_boxplots_robustSimStudy

% performance
s5 = subplot('Position',[0.4,0.08,0.5,0.15]);
boxplot(converged./t_cpus*60,'OutlierSize',4,'Symbol','k.','Width',0.8);hold on;
s5.YScale = 'log';
s5.TickDir = 'out';
s5.XTick = [];
s5.FontSize = 8;
if labelFlag
    ylabel('performance')
else
    s5.YTickLabel = '';
end
s5.YTick = [1e-1,1e0,1e1,1e2];
ylim([5e-2,1e2])

color_boxplots_robustSimStudy
