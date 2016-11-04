%% figure5
orig_bp=pwd;
example5
bp = '../..';
cd(bp);


%% this part shows the dependence on the ratio of the mus
load('contract_meas')

ax1 = subplot('Position',[0.565,0.915-exw*Dy/Dx_-0.225,0.39,0.215]);
unis = unique(-allp(:,7)./allp(:,3));
cols = parula(length(unis)+1);
for muc = unique(-allp(:,3))'
    for ups = unique(allp(:,7))'
        subind = allp(:,3)==-muc&allp(:,7)==ups&allp(:,6)==10;
        coli = find(unis==ups/muc);
        x = -allp(subind,end);
        if(length(x)>2)
            w0 = repmat(allw(subind,1),1,size(allw(subind,:),2));
            y = max(-(allw(subind,:)-w0)./w0,[],2);
            [sortedX, sortIndex] = sort(x);
            plot(sortedX,y(sortIndex),'DisplayName',['\upsilon/\mu = ' num2str(ups/muc)],'Color',cols(coli,:))
            hold on
        end
    end
end
colormap(ax1,cols)
cb = colorbar('Location','south','Ticks',linspace(0,1,5),'TickLabels',{'0.3','1.0','3.3','10','33'});
title(cb,'\upsilon/\mu_c');
pos = cb.Position;
cb.Position = [pos(1)+pos(3)/2 pos(2) pos(3)/2 pos(4)/3];
xlabel('Stiffness Asymmetry (\mu_e/\mu_c)')
ylabel('Maximum Strain')

annotation('textbox', [0.01 0.89 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.01 0.68 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.49 0.68 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

cd('../figures')
print('-depsc','-r0',['figure5.eps']);
cd(orig_bp)