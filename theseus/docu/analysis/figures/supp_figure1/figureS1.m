orig_bp=pwd;
bp = '../../../data/';
cd(bp);
load('extend_meas2')
h2=figure;


ax1 = subplot('Position',[0.1,0.1,0.35,0.3]);
ax2 = subplot('Position',[0.55,0.1,0.35,0.3]);
for ind=1:size(allt,1)
    
    r = allp(ind,10)*10;
    tr = 1/r;
    L = allp(ind,2);
    lc = allp(ind,5);
    mu = abs(allp(ind,3));
    xi = allp(ind,6)/10;
    ups = allp(ind,7);
    
    sig = abs(allp(ind,11));
    
    G_pred = 2*mu/lc;
    
    t = allt(ind,:);
    tstop = find(t==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(t);
    end
    
    myt = t(1:tstop)/10;
    myg = allg(ind,1:tstop)*10;
    sl2 = cumtrapz(myt,myg,2);
    sl = diff(log(sl2),1,2)./diff(log(myt),1,2);

    cutoff = find(sl==min(sl),1);
    term = cutoff+find(sl2(cutoff:end)>0.15,1);
    if(isempty(term))
        term = length(sl);
    end
    spt = cutoff+find(sl(cutoff:term-1)>0.25,1);
    spt
    if(~isempty(spt)&&allp(ind,6)==1)
        G = abs(allp(ind,11))./sl2(spt);
        if(myt(spt)>50)
        axes(ax1);
        plot(myt(1:spt),sl2(1:spt)/sig*G_pred)
        hold on
        end
        axes(ax2);
        plot(G_pred,G,'.','Color',[0.25 0.25 0.25])
        hold on
        
    end
    
end
axes(ax2);
xlabel('Estimated elastic modulus (2\mu/l_c)')
ylabel('Measured elastic modulus, G_0 (nN/\mum)')
axes(ax1);
ylabel('Normalized Strain (\gamma/\sigma\cdot2\mu/l_c)')
xlabel('Time (s)')

annotation('textbox', [0.005 0.37 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.46 0.37 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

cd('../figures')
print('-depsc','-r0','figureS1.eps');
cd(orig_bp)
