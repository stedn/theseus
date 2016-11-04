

%% figureS3
orig_bp = pwd;

exampleS3_1
exampleS3_2



annotation('textbox', [0.01 0.89 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.01 0.68 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

annotation('textbox', [0.00 0.45 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.48 0.45 0.05 0.05],'String','d)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
cd('../../../figures')
print('-depsc','-r0',['figureS3.eps']);
cd(orig_bp)





