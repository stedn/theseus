bp = '../../../data/';
cd(bp);

h2=figure;


load('domainllc_meas')


subplot('Position',[0.075 0.95-0.4 0.35 0.25])


for ind=1:size(allp,1)
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10
    ups = allp(ind,7)
    tscale=L*xi/(ups*mu)^(1/2);
    sscale=(ups*mu)^(1/2)/lc;
    nscale = (L/lc -1)^2*xi;

    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop)/10;
    g = allg(ind,floor(tstop/2):tstop);

    if(t(end)>tscale*10)
         plot(L,mean(g)*xi,'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10) '  \xi = ' num2str(allp(ind,6)) '  \upsilon = ' num2str(allp(ind,7)) '  L = ' num2str(allp(ind,2)) '  l_c = ' num2str(allp(ind,5)) '  t_{last} = ' num2str(t(end))])
         hold on
    end
end
xlabel('Filament Length , L (\mum)')
ylabel('Normalized Strain Rate ($$\dot{\gamma}\xi$$)','interpreter','latex')
ylim ([0 0.004])





subplot('Position',[0.575 0.95-0.4 0.35 0.25])
for ind=1:size(allp,1)
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10
    ups = allp(ind,7)
    tscale=L*xi/(ups*mu)^(1/2);
    sscale=(ups*mu)^(1/2)/lc;
    nscale = (L/lc -1)^2*xi;

    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    g = allg(ind,floor(tstop/2):tstop);

    if(t(end)>tscale*10)
         plot(lc,mean(g)/ups*xi*L^1.1,'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10) '  \xi = ' num2str(allp(ind,6)) '  \upsilon = ' num2str(allp(ind,7)) '  L = ' num2str(allp(ind,2)) '  l_c = ' num2str(allp(ind,5)) '  t_{last} = ' num2str(t(end))])
         hold on
    end
end
xlabel('Cross-link spacing, l_c (\mum)')
ylabel('Normalized Strain Rate ($$\dot{\gamma}\xi L^{1.1} $$)','interpreter','latex')
ylim ([0 .11])

cd('../figures')
print('-depsc','-r0',['figure10.eps']);
