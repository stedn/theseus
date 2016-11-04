orig_bp=pwd
bp = '../../../data/';
cd(bp);

%% load premeasured data
load('extend_meas2')
h2=figure;


%% subplot for effective viscosity scaling

subplot('Position',[0.1,0.1,0.35,0.3]);
for ind=1:size(allt,1)

    t = allt(ind,:)/10;

    L = allp(ind,2)
    lc = allp(ind,5)
    mu = abs(allp(ind,3))
    xi = allp(ind,6)/10

    tstop = find(t==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(t);
    end

    myt = t(1:tstop);
    myg = allg(ind,1:tstop);
    sl2 = cumtrapz(myt,myg,2);
    sl = diff(log(sl2),1,2)./diff(log(myt),1,2);

    cutoff = find(sl==min(sl),1);
    term = cutoff+find(sl2(cutoff:end)>0.15,1);
    if(isempty(term))
        term = length(sl);
    end
    spt = cutoff+find(sl(cutoff:term-1)>0.75,1);
    if(~isempty(spt))
        eta = abs(allp(ind,11))./mean(allg(ind,spt:end),2);
        plot(L/lc-1,eta/xi,'.','Color',[0.25 0.25 0.25])
        hold on

    end

end
xlabel('Number of crosslinks per filament (L/l_c-1)')
ylabel('Normalized effective viscosity (\eta/\xi)')
plot(linspace(0,25,20),3*pi*linspace(0,25,20).^2,'k:')










subplot('Position',[0.55,0.1,0.35,0.3]);
for ind=1:size(allt,1)

    t = allt(ind,:)/10;

    L = allp(ind,2)
    lc = allp(ind,5)
    mu = abs(allp(ind,3))
    xi = allp(ind,6)/10

    tstop = find(t==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(t);
    end

    myt = t(1:tstop);
    myg = allg(ind,1:tstop);
    sl2 = cumtrapz(myt,myg,2);
    sl = diff(log(sl2),1,2)./diff(log(myt),1,2);

    cutoff = find(sl==min(sl),1);
    term = cutoff+find(sl2(cutoff:end)>0.15,1);
    if(isempty(term))
        term = length(sl);
    end
    spt = cutoff+find(sl(cutoff:term-1)>0.75,1);
    if(~isempty(spt))
        tau=myt(cutoff+find(sl(cutoff:term-1)>0.8,1));
        if(~isempty(tau)&&tau<t(end))
            loglog(L^2/lc/mu.*xi,tau,'.','Color',[0.25 0.25 0.25])
            hold on
        end
    end

end
loglog([10 1000],[10 1000],'k:')
xlim([10 1000])
ylim([10 1000])
xlabel('Predicted \tau_c (L^2\xi/l_c\mu_e)')
ylabel('Simulated Transition Time \tau_c (s)')








annotation('textbox', [0.005 0.37 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.46 0.37 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

cd('../figures')
print('-depsc','-r0','figure3.eps');
cd(orig_bp)