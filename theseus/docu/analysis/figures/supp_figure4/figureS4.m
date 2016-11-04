exampleS4_1
exampleS4_2

bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);

annotation('textbox', [0.005 0.92 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.5 0.92 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

print('-depsc','-r0',['figureS4.eps']);
