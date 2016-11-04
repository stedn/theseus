%% figure2S
h2=figure;
orig_bp=pwd;
bp = '../../../data/';
cd(bp);
load('contract_meas')




subplot('Position',[0.075 0.95-0.2 0.35 0.2])
for ind=1:size(allt,1)
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10
    ups = allp(ind,7)
    
    t = allt(ind,:)/10;
    w = allw(ind,:);
    w = (w(1)-w)/w(1);
    taui = find(w>0.7*max(w));
    stotau = t(taui(1));
    loglog(L*xi/ups,stotau,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6))])
    hold on
end
loglog([0.2,200],[0.2,200],'k:')
xlim([0.2,200])
ylim([0.2,200])

xlabel('Predicted \tau_m = (L\xi/\upsilon)')
ylabel('Time of Max Strain')


load('astress_meas')

subplot('Position',[0.55 0.95-0.2 0.35 0.2])
for ind=1:size(allt,1)
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10
    ups = allp(ind,7)
    phi = allp(ind,8)
    
    mx = max(allf(ind,:));
    if(allp(ind,end)>0)
        plot(phi,mx*lc,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6)) ' ' num2str(allp(ind,end))])
        hold on
    end
end

xlabel('Activity Fraction $$\phi$$','interpreter','latex')
ylabel('Normalized Max Stress ($$\sigma \cdot l_c$$)','interpreter','latex')



annotation('textbox', [0.02 0.91 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.42 0.91 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

cd('../figures')
print('-depsc','-r0',['figureS2.eps']);
cd(orig_bp)