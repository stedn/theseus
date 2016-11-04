orig_bp=pwd;
example6

bp = '../..';
cd(bp);
load('astress_meas_correct')

ax2 = subplot('Position',[0.565,0.915-exw*Dy/Dx_-0.225,0.39,0.215]);
for ind=1:size(allt,1)
    t = allt(ind,:)/10;

    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10

    mx = max(allf(ind,:));
    if(mx>0.005&&allp(ind,end)>0)
        taui = find(allf(ind,:)>0.9*mx);
        stotau = t(taui(1));
        tau_pred = L*xi/sqrt(ups*mu)
        loglog(tau_pred,stotau,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6)) ' ' num2str(allp(ind,end))])
        hold on
    end
end
loglog([0.1,1000],[0.1,1000],'k:')
% xlim([0.1 1000])
% ylim([0.1 1000])
xlabel('$$\tau_a = $$ ($$L\xi/\sqrt{\mu_e\upsilon}$$)','interpreter','latex')
ylabel('Time of Max Stress')

annotation('textbox', [0.01 0.89 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.00 0.68 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.49 0.68 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

cd('../figures')
print('-depsc','-r0',['figure6.eps']);
cd(orig_bp)