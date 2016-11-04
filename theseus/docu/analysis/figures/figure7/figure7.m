orig_bp=pwd;
h2=figure;
example7_1
example7_2
example7_3
example7_4

bp = '../..';
cd(bp);
load('actrec_meas')

subplot('Position',[0.525 topp-0.3 0.4 0.275])

indabl = find(allp(:,6)==10&allp(:,5)==0.3&allp(:,7)==0.1&allp(:,8)==0.25&...
    (allp(:,10)==0.0001|allp(:,10)==0.001|allp(:,10)==0.01|allp(:,10)==0.1));
[dum,srt] = sort(allp(indabl,10));
for ind=indabl(srt)'

    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = abs(allp(ind,3))
    xi = allp(ind,6)/10

    alln{ind}
    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop)/10;
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);

    if(1)
                plot([0 t],[0 sl],'DisplayName',['\tau_r = ' num2str(tr,4)])
                hold on
    end
end
ylabel('Stress (nN)')
xlabel('Time (s)')
legend('Location','northeast')
xlim([0 250])

indabl = find(allp(:,6)==10&allp(:,5)==0.3&allp(:,7)==0.1&allp(:,8)==0.25);
    %(allp(:,10)==0.0003|allp(:,10)==0.003|allp(:,10)==0.03|allp(:,10)==0.3|allp(:,10)==3));
[dum,srt] = sort(allp(indabl,10));
subplot('Position',[0.1 topp-0.3*2 0.3 0.2])
st_x=[];
st_y=[];
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

    st_x=[st_x tr/tscale];
    st_y=[st_y mean(sl(end-10:end))];
end
semilogx(st_x,st_y,'Color',[0.25,0.25,0.25])

ylabel('Steady State Stress (nN)')
xlabel('Normalized Recycling Time (\tau_r/\tau_a)')
xlim([0.0001 100])
set(gca,'XTick',[0.01 1 100],'XTickLabel',[0.01 1 100])




load('actrec_meas')


subplot('Position',[0.525+0.025 topp-0.3*2 0.4-0.05 0.23])
for ind=1:size(allt,1)
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10
    ups = allp(ind,7)

    tscale=L*xi/(ups*mu)^(1/2);
    sscale=0.1*(ups*mu)^(1/2)/lc;

    tstop = find(allt(ind,:)==0,2);
    if(length(tstop)>1)
        tstop = tstop(2)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop)/10;
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);
    strt = find(abs(t-2*tr)==min(abs(t-2*tr)))
    strss = mean(sl(strt:end))
    loglog(tr/tscale,strss/sscale,'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau = ' num2str(tr) ',  \xi = ' num2str(xi) ',  \upsilon = ' num2str(ups)])
    hold on

end
myx = logspace(-4,4,35);
loglog(myx,1./(1./myx+myx),'--')
ylabel('Normalized Steady State Stress (\sigma/\sigma_a)')
xlabel('Normalized Recycling Time (\tau_r/\tau_a)')

annotation('textbox', [0.005 0.66 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.44 0.66 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.001 0.27 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.44 0.3 0.05 0.05],'String','d)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

cd('../figures')
print('-depsc','-r0',['figure7.eps']);
cd(orig_bp)