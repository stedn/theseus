orig_bp=pwd;
example4_1
example4_2



%% plot the strain curves for the recycling and the no recycling case
subplot('Position',[0.55 0.94-0.2*Dy/Dx_*2 0.4 0.2*Dy/Dx_*2])


hold on
ylabel('Strain')
xlabel('Time (s)')
bp = '../..';
cd(bp);
load('extendrec_meas')

indabl = find(allp(:,6)==1&allp(:,5)==0.5&allp(:,11)==-0.0002&allp(:,12)==75);
[dum,srt] = sort(allp(indabl,10));
for ind=indabl(srt)'

    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop)/10;
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop)*10;

    r = allp(ind,10)*10

    if(1)
                plot(t,cumtrapz(t,g,2),'DisplayName',['\tau_r = ' num2str(1/r) ])
     hold on
    end
end


%% just fyi this for loop is only going to plot one thing

load('extend_meas')
indabl = find(allp(:,6)==1&allp(:,5)==0.5&allp(:,11)==-0.0002&allp(:,2)==5);
[dum,srt] = sort(allp(indabl,10));
for ind=indabl(srt)'

    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop)/10;
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop)*10;

    if(1)
                plot(t,cumtrapz(t,g,2),'DisplayName','\tau_r = \infty' )
     hold on
    end
end
xlim([0 200])
legend('Location','northeast')





%% this is going to plot all of the strain rates for the experiments above


load('extendrec_meas')
subplot('Position',[0.075+0.0125 0.94-0.25*2 0.375 0.2])
st_x=[];
st_y=[];

indabl = find(allp(:,6)==1&allp(:,5)==0.5&allp(:,11)==-0.0002&allp(:,2)==5);
[dum,srt] = sort(allp(indabl,10));

for ind=indabl(srt)'
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = abs(allp(ind,3))
    xi = allp(ind,6)/10

    sig = abs(allp(ind,11))

    ntscale = L^2*xi/lc/mu;

    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end

    t = allt(ind,1:tstop)/10;
    g = allg(ind,floor(tstop/2):tstop);


    st_x=[st_x tr/ntscale];
    st_y=[st_y sig/mean(g)];

end

load('extend_meas')

indabl = find(allp(:,6)==1&allp(:,5)==0.5&allp(:,11)==-0.0002&allp(:,2)==5);
[dum,srt] = sort(allp(indabl,10));

for ind=indabl(srt)'
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = abs(allp(ind,3))
    xi = allp(ind,6)/10

    sig = abs(allp(ind,11))

    ntscale = L^2*xi/lc/mu;

    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end

    t = allt(ind,1:tstop)/10;
    g = allg(ind,floor(tstop/2):tstop);

    st_x=[10^4.7 st_x ];
    st_y=[sig/mean(g) st_y ];
end

semilogx(st_x,st_y,'Color',[0.25,0.25,0.25])
hold on
semilogx(st_x(1),st_y(1),'o','Color',[0.25,0.25,0.25],'MarkerSize',6)
xlim([0.000005,10^5])
ylabel('Effective Viscosity (nNs/\mum)')
xlabel('Normalized Recycling Time (\tau_r/\tau_c)')



%% and now the big one where we plot everything under the sun scaled

load('extendrec_meas')

subplot('Position',[0.55+0.0125 0.94-0.25*2 0.375 0.2])
bigtrack=[]
for ind=1:size(allt,1)
    r = allp(ind,10)*10;
    tr = 1/r;
    L = allp(ind,2);
    lc = allp(ind,5);
    mu = 100*abs(allp(ind,3));
    xi = allp(ind,6)/10;

    sig = abs(allp(ind,11));

    ntscale = L^2*xi/lc/mu*4*pi;
    nscale = (L/lc-1)^2*xi*4*pi;

    tstop = find(allt(ind,:)==0,2);
    if(length(tstop)>1)
        tstop = tstop(2)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop)/10;
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);

    if(allp(ind,11)<0&&tr>=xi&&sig>1e-5)
         gam = mean(g(find(abs(t-2*tr)==min(abs(t-2*tr))):end))
                loglog(tr/ntscale,sig/gam/nscale,'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau = ' num2str(tr) ',  \xi = ' num2str(xi) ',  \sigma = ' num2str(sig)])
     hold on
        if(sig/gam>nscale*5)
            bigtrack
        end
    end
end


%% and this just adds the no recycling data
load('extend_meas2')

for ind=1:size(allt,1)
    r = allp(ind,10)*10;
    tr = 1/r;
    L = allp(ind,2);
    lc = allp(ind,5);
    mu = 100*abs(allp(ind,3));
    xi = allp(ind,6)/10;

    sig = abs(allp(ind,11));

    ntscale = L^2*xi/lc/mu*4*pi;
    nscale = (L/lc-1)^2*xi*4*pi;

    tstop = find(allt(ind,:)==0,2);
    if(length(tstop)>1)
        tstop = tstop(2)-1;
    else
        tstop = length(allt(ind,:));
    end

    t = allt(ind,1:tstop)/10;
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
    if(~isempty(spt))
        eta=sig/mean(allg(ind,spt:end),2);
        loglog(10^4.7+100*randn,eta/nscale,'o','Color',[0.25 0.25 0.25],'MarkerSize',6,'DisplayName',['\tau = ' num2str(1/allp(ind,10)) ',  \xi = ' num2str(allp(ind,6)) ',  \upsilon = ' num2str(allp(ind,7))])
        hold on
    end
end
myx = logspace(-5,5,35);
loglog(myx,1./(1+1./myx.^(3/4)),'--')

loglog([0.00001 100000],[1 1],':','Color',[0.25 0.25 0.25])
ylim([0.001 8])
ylabel('Normalized Viscosity (\eta/\eta_c)')
xlabel('Normalized Recycling Time (\tau_r/\tau_c)')



%% drop some nnotations
annotation('textbox', [0.00 0.915 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.49 0.915 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.00 0.61 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.47 0.61 0.05 0.05],'String','d)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])


cd('../figures')
print('-depsc','-r0',['figure4.eps']);
cd(orig_bp)