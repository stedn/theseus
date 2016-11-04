
bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);

load('domainllc_meas')


indabl = find(allp(:,6)>0);
[dum,srt] = sort(allp(indabl,10));

for ind=indabl(srt)'
    mu = 100*abs(allp(ind,3));
%     if(allp(ind,3)<0)
%         mu=100*mu;
%     end
    tr = 1./allp(ind,10);
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
    sscale=10*(allp(ind,7)*mu)^(1/2)/allp(ind,5);
    nscale = (allp(ind,2)/allp(ind,5))^2*allp(ind,6);
    ntscale = allp(ind,2).^2.*allp(ind,6)./allp(ind,5)./mu;

    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    g = allg(ind,floor(tstop/2):tstop);

    if(t(end)>0)
         semilogx(tr/tscale,mean(g)/allp(ind,7)*allp(ind,6)*allp(ind,2),'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10) '  \xi = ' num2str(allp(ind,6)) '  \upsilon = ' num2str(allp(ind,7)) '  t_{last} = ' num2str(t(end))])
         hold on
    end
end
xlabel('Normalized Recycling Time (\tau_r/\tau_a)')
ylabel('Normalized Strain Rate ($$\dot{\gamma}\xi/\upsilon$$)','interpreter','latex')

print('-depsc','-r0',['figure6c.eps']);



load('domain_meas')

h2=figure
for ind=1:size(allp,1)
    mu = 100*abs(allp(ind,3));
    alln{ind}
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
    ntscale = allp(ind,2)^2*allp(ind,6)/allp(ind,5)/abs(allp(ind,3));
    sscale=1/10*(allp(ind,7)*mu)^(1/2)/allp(ind,5);
    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);

    if(1)
                semilogx(tr/tscale,mean(sl(find(abs(t-2*tr)==min(abs(t-2*tr))):end))/sscale,'.','DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10,3)])
                hold on
    end
end
