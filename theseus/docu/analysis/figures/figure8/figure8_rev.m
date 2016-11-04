orig_bp=pwd;
example8_1
example8_2

bp = '../..';
cd(bp);
load('domain_meas')

subplot('Position',[exw+0.15 0.95-(0.8-exw)*Dy/Dx_*1.525 0.8-exw (0.8-exw)*Dy/Dx_*1.5])
indabl = find(allp(:,6)==10&allp(:,7)==0.1);
[dum,srt] = sort(allp(indabl,10));
for ind=indabl(srt)'
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10
    ups = allp(ind,7)

    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop)/10;
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);

    if(t(end)>200)
                plot(t,cumtrapz(t,g,2),'DisplayName',['\tau_r = ' num2str(tr,4)])
     hold on
    end
end
legend('Location','northwest')
xlabel('Time (s)')
ylabel('Strain')


load('domain_meas')


indabl = find(allp(:,6)==10&allp(:,7)==0.1);
[dum,srt] = sort(allp(indabl,10));
subplot('Position',[0.1 0.9-0.3*2 0.3 0.22])
example8_1_spatial
example8_2_spatial

subplot('Position',[exw+0.2 0.9-0.3*2 0.8-exw-0.07 0.25])
st_x = []
st_y = []
for ind=indabl(srt)'
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10
    ups = allp(ind,7)
    tscale=L*xi/(ups*mu)^(1/2);

    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop)/10;
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);

    if(t(end)>2*tscale)
                st_x=[st_x tr/tscale];
                st_y=[st_y mean(g(end-100:end))];

    end
end
semilogx(st_x,st_y,'Color',[0.25,0.25,0.25])
xlabel('Normalized Recycling Time (\tau_r/\tau_a)')
ylabel('Strain Rate (1/s)','interpreter','latex')





annotation('textbox', [0.005 0.91 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.44 0.91 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.005 0.9-0.4-0.01 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.42 0.9-0.4+0.01 0.05 0.05],'String','d)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

cd('../../../figures')
print('-depsc','-r0',['figure8.eps']);
cd(orig_bp)