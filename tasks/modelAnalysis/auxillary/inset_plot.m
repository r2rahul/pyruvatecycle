function inset_plot(lgbar,lgerbar,nadlgbar,nadlgerbar,fname)
%cmap=[0.6,0.55686,0.76471;0.9451,0.63922,0.25098];
%cmap=[0.1451,0.1451,0.1451;0.3882,0.3882,0.3882];
%cmap=[0.1451,0.1451,0.1451;0.58824,0.58824,0.58824];
c1 = [37,37,37] ./ 255;
c2 = [99,99,99] ./ 255;
cmap = [c1; c2];
%seq1=[0.74118,0.74118,0.74118];
seq1 = [217,217,217] ./ 255;
h1=figure(1); clf;
colormap(cmap)
bb=bar(lgbar);
set(bb,'BarWidth',1,'EdgeColor','none'); % The bars will now touch each other
set(gca,'FontSize',12,'TickLength',[0.0,0.05],'LineWidth',0.2,'box','off'...
    ,'TickDir','out');

%title('3mM Glucose - siControl Experimental and Simulation Values');
%xlabel('Metabolite Name')
set(gca,'XLim',[0 8],'XTick',1:7);
set(gca,'XTickLabel',{'Pyruvate','Citrate','AKG','Succinate','Fumarate','Malate','Lactate'})
xtickangle(45)
%set(gca,'Xcolor',[0.7412,0.7412,0.7412]);
%set(gca,'Ycolor',[0.7412,0.7412,0.7412]);
ylabel('Metabolite Concentration in (M)')
%legend('Simulation','Observed',2,'Location','NorthEast')
%gridxy([],get(gca,'ytick'),'color',[0.7412,0.7412,0.7412],'linewidth',1)
%set(gca,'Layer','bottom')
%drawGrid(gca,0.1)
hold on;
numgroups =size(lgbar, 1);
numbars =size(lgbar, 2);
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, lgbar(:,i), lgerbar(:,i), 'k', 'linestyle', 'none');
end
set(gca,'YGrid','on')
set(gca,'GridLineStyle','-')
set(gca,'XColor',seq1);
set(gca,'YColor',seq1);
Caxes = copyobj(gca,gcf);
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','on');
%hold on;
box off;
myAxes = axes('parent', h1, 'Position', [0.45, 0.45, 0.2, 0.2]);
bb=bar(nadlgbar);
set(bb,'BarWidth',1,'EdgeColor','none'); % The bars will now touch each other
set(gca,'FontSize',12,'box','off');
set(gca,'XLim',[0 2],'XTick',1,'TickLength',[0.0,0.025]);
set(gca,'XTickLabel',{'NADPH'})
xtickangle(45)
hold on
numgroups = size(nadlgbar, 1);
numbars = size(nadlgbar, 2);
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, nadlgbar(:,i), nadlgerbar(:,i), 'k', 'linestyle', 'none');
end
set(gca,'YGrid','on')
set(gca,'GridLineStyle','-')
%[0.7412,0.7412,0.7412]
set(gca,'XColor',seq1);
%[0.7412,0.7412,0.7412]
set(gca,'YColor',seq1);
Caxes = copyobj(gca,gcf);
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','on');
lgd = legend('Simulation', 'Observed', 'Location', 'best', 'NumColumns', 1);
legend boxoff
lgd.Position =[0.45, 0.75, 0.2, 0.2];
%applyhatch(h1,'\-x.+');
print(h1,'-dsvg','InsetFigures');
print(h1,'-dtiff',fname);
print(h1,'-dpng',fname);
hold off