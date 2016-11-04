figure
orig_bp=pwd;
bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('extendrec_meas')

subplot('Position',[0.125 0.95-0.2 0.325 0.2])
for ind=1:size(allt,1)
    
       
    r = allp(ind,10)*10;
    tr = 1/r;
    L = allp(ind,2);
    lc = allp(ind,5);
    mu = 100*abs(allp(ind,3));
    xi = allp(ind,6)/10;
    ups = allp(ind,7);
    
    ntscale=L*xi/ups;
    sscale=(ups*mu)^(1/2)/lc;
    nscale = (L/lc -1)^2*xi;
    
    tstop = find(allt(ind,:)==0,2);
    if(length(tstop)>1)
        tstop = tstop(2)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);
    
    if(allp(ind,11)<0&&tr>=allp(ind,6))
                loglog(tr/ntscale,sig/mean(g(find(abs(t-2*tr)==min(abs(t-2*tr))):end))/nscale,'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau = ' num2str(1/allp(ind,10)) ',  \xi = ' num2str(allp(ind,6)) ',  \sigma = ' num2str(allp(ind,11))])
     hold on
    end
end

load('extend_meas2')

subplot('Position',[0.125 0.95-0.2 0.325 0.2])
for ind=1:size(allt,1)
    mu = abs(allp(ind,3));
    
    nscale = pi/4*(allp(ind,2)/allp(ind,5)-1)^2*allp(ind,6);
    ntscale = allp(ind,2)^2*allp(ind,6)/allp(ind,5)/abs(allp(ind,3));
    
    tstop = find(allt(ind,:)==0,2);
    if(length(tstop)>1)
        tstop = tstop(2)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);
    
    sl2 = cumtrapz(t,g,2);
    sl = diff(log(sl2),1,2)./diff(log(t),1,2);

    cutoff = find(sl==min(sl),1);
    term = cutoff+find(sl2(cutoff:end)>0.15,1);
    if(isempty(term))
        term = length(sl);
    end
    spt = cutoff+find(sl(cutoff:term-1)>0.75,1);
    eta=abs(allp(ind,11))./mean(allg(ind,spt:end),2);
    if(~isempty(spt))
                loglog(10^4.7+100*randn,eta/nscale,'o','Color',[0.25 0.25 0.25],'MarkerSize',6,'DisplayName',['\tau = ' num2str(1/allp(ind,10)) ',  \xi = ' num2str(allp(ind,6)) ',  \upsilon = ' num2str(allp(ind,7))])
     hold on
    end
end
loglog([0.00001 100000],[1 1],':','Color',[0.25 0.25 0.25])
ylim([0.001 5])
ylabel('Normalized Viscosity (\eta/\eta_c)')
xlabel('Normalized Recycling Time (\tau_r/\tau_c)')

bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('actrec_meas')


subplot('Position',[0.575 0.95-0.2 0.325 0.2])
for ind=1:size(allt,1)
    mu = 100*abs(allp(ind,3));
    
    tr = 1./allp(ind,10);
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
    tscale2=allp(ind,2)*allp(ind,6)/allp(ind,7);
    sscale=1/10*(allp(ind,7)*mu)^(1/2)/allp(ind,5);
    
    tstop = find(allt(ind,:)==0,2);
    if(length(tstop)>1)
        tstop = tstop(2)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);
    
    if(1)
                loglog(tr/tscale,mean(sl(find(abs(t-2*tr)==min(abs(t-2*tr))):end))/sscale,'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau = ' num2str(1/allp(ind,10)) ',  \xi = ' num2str(allp(ind,6)) ',  \upsilon = ' num2str(allp(ind,7))])
     hold on
    end
end
ylabel('Normalized Steady State Stress (\sigma/\sigma_a)')
xlabel('Normalized Recycling Time (\tau_r/\tau_a)')

annotation('textbox', [0.01 0.91 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.46 0.91 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

cd('../figures')
print('-depsc','-r0',['figure4S.eps']);
cd(orig_bp)