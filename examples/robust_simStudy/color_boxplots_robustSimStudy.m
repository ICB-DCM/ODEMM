box off
boxes = findobj(gca,'Tag','Box');
median = findobj(gca,'Tag','Median');
patch(get(boxes(4),'XData'),get(boxes(4),'YData'),color.norm,'FaceAlpha',1); hold on;
patch(get(boxes(3),'XData'),get(boxes(3),'YData'),color.skew_norm,'FaceAlpha',1); hold on;
patch(get(boxes(2),'XData'),get(boxes(2),'YData'),color.students_t,'FaceAlpha',1); hold on;
patch(get(boxes(1),'XData'),get(boxes(1),'YData'),color.neg_binomial,'FaceAlpha',1); hold on;
plot(get(median(1),'XData'),get(median(1),'YData'),'k-','Linewidth',0.5); hold on;
plot(get(median(2),'XData'),get(median(2),'YData'),'k-','Linewidth',0.5);
plot(get(median(3),'XData'),get(median(3),'YData'),'k-','Linewidth',0.5);
plot(get(median(4),'XData'),get(median(4),'YData'),'k-','Linewidth',0.5);

tmpline = findobj(gca,'Type','line');
for indLine = length(tmpline):-1:(length(tmpline)-8)
    tmpline(indLine).LineStyle = '-';
end