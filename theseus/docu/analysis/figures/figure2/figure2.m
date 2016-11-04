orig_bp = pwd
bp = '../../../data/examples/fig3/';
code = 'zkyqjony';%gcqbbcyr
cd(bp)

makemovs = 0;

%% load param file and decipher params
fid = fopen([bp code '_scr.txt']);
C = textscan(fid, '%s','delimiter', '\n');
fclose(fid);
pare = strsplit(C{1}{9}, '>');
paree = strsplit(pare{1}, ' ');
paree = {paree{2:end}};
zet=str2num(paree{2});L=str2num(paree{3});mu=-str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6});
xi=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
r=str2num(paree{11});sig=str2num(paree{12});Dx=str2num(paree{13});Dy=str2num(paree{14});Df=str2num(paree{15});
Dw=str2num(paree{16});ls=str2num(paree{17});lf=str2num(paree{18});

xi=xi/10
r = r*10

%% load simulation data
A = importdata([bp code '_out.txt']);
A = A.data;
if(size(A,1)==1)
    imp2 = importdata([bp code '_out.txt'],' ',9);
    if(isfield(imp2,'data'))
        A = [A;imp2.data];
    end
end
t = A(:,1)/10;
zt = A(:,2:end);

%% store initial positions and initial measurements (all 0)
op = reshape(zt(1,:),[],2);
tl=0;
stof = 0;
stog = 0;
stoa = 0;
stot = 0;


%% setup timepoints and space points to measure
inds = 1:ceil(size(zt,1)/250):size(zt,1);
inds = inds(2:end);
ex_indis = [1 40 200];
bpos = linspace(0,Dx,51);
bpos = bpos(1:end-1)+bpos(2)/2;
ll = 4;
rl = 16;

h2 = figure;
if(makemovs)
    clear mov mov2
    h1 = figure;
    h = figure('Position', [50, 100, 100+600*Dx/Dy, 600]);
end

%% just setup stupid plotting colors
lst = size(zt,1);
trp = repmat((1:lst)'/lst,1,3);
temp=flipud(winter(lst));
temp2 = bone(2*lst);
cc = (1-trp.^2).*temp2(lst+1:end,:)+(trp.^2).*temp(1:lst,:);
temp=hot(2*lst);
cc2 = (1-trp.^2).*temp2(lst+1:end,:)+(trp.^2).*temp(1:lst,:);
indi = 1;


Dx_ = 0.5*Dx;
Dx__ = 0.35*Dx;

%% color limits for strains
cdmn=0.06;
edmn=0.06;

%% figure width
www = 0.37;


%% for loop over timepoints to display
for ind = inds
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),Dx),mod(p(:,2),Dy)];
    if(makemovs)
        figure(h)
        clf
        netplot_str(p,L,lf,ls,Dx,Dy,cc,cc2,edmn,cdmn);
        xlim([0 Dx_])
        set(gca,'xtick',[],'ytick',[],'box','on')

        colormap([flipud(cc2);cc])
        colorbar('Location','westoutside','box','on','Ticks',[0 0.5 1],'TickLabels',{num2str(-cdmn), '0.00', num2str(edmn)},'TickDirection','out');
        drawnow
        mov(indi) = getframe(h);
    end
    figure(h2)
    if(indi==ex_indis(1))
        subplot('Position',[0.06 0.97-www*Dy/Dx_ www www*Dy/Dx_*0.9])
        pat=patch([Dx__-Dx*Dw Dx__-Dx*Dw Dx__ Dx__],[0 Dy Dy 0],[.7 .5 0]);
        set(pat,'FaceAlpha',0.25,'EdgeColor','none');
        netplot_str(p,L,lf,ls,Dx,Dy,cc,cc2,edmn,cdmn);
        ylabel(['t = ' num2str(t(ind)) ' s'])
        xlim([0 Dx_])
        set(gca,'xtick',[],'ytick',[],'box','on')

        colormap([flipud(cc2);cc])
        cb=colorbar('Location','east','box','on','Ticks',[0 0.5 1],'TickLabels',{num2str(-cdmn), '0.00', num2str(edmn)},'TickDirection','out');
        ylabel(cb,'Strain')
        pos = cb.Position;
        cb.Position = [pos(1)+pos(3)/8 pos(2)+pos(4)/2 pos(3)/2 pos(4)/2];
    elseif(indi==ex_indis(2))
        subplot('Position',[0.06 0.97-www*Dy/Dx_*2 www www*Dy/Dx_*0.9])
        pat=patch([Dx__-Dx*Dw Dx__-Dx*Dw Dx__ Dx__],[0 Dy Dy 0],[.7 .5 0]);
        set(pat,'FaceAlpha',0.25,'EdgeColor','none');
        netplot_str(p,L,lf,ls,Dx,Dy,cc,cc2,edmn,cdmn);
        ylabel(['t = ' num2str(t(ind)) ' s'])
        xlim([0 Dx_])
        set(gca,'xtick',[],'ytick',[],'box','on')

    elseif(indi==ex_indis(3))
        subplot('Position',[0.06 0.97-www*Dy/Dx_*3 www www*Dy/Dx_*0.9])
        pat=patch([Dx__-Dx*Dw Dx__-Dx*Dw Dx__ Dx__],[0 Dy Dy 0],[.7 .5 0]);
        set(pat,'FaceAlpha',0.25,'EdgeColor','none');
        netplot_str(p,L,lf,ls,Dx,Dy,cc,cc2,edmn,cdmn);
        ylabel(['t = ' num2str(t(ind)) ' s'])
        xlim([0 Dx_])
        set(gca,'ytick',[],'xtick',[],'box','on')



    end

    dp = (p-op);
    op = p;

    % remove data if has moved farther than realistically possible
    % these events are due to crossing domain boundary or recycling
    jumpcut = 50*median(abs(dp(1:1:end,1)));
    subind = abs(dp(1:2:end-1,2))>jumpcut|abs(dp(2:2:end,2))>jumpcut|abs(dp(1:2:end-1,1))>jumpcut|abs(dp(2:2:end,1))>jumpcut;
    subind = [subind subind]';
    subind = subind(:);
    dp(subind,:)=[];
    p(subind,:)=[];

    % compute velocities
    v = dp/(t(ind)-tl);
    tl = t(ind);

    % compute filament strain
    [XY,sx,sy]=get_str(p,L,lf,ls,Dx,Dy);

    % compute filament tension
    fx = mu*sx;
    fy = mu*sy;
    if(mu<0)
        fx = -fx.*(1+99*double(sx>0));
        fy = -fy.*(1+99*double(sy>0));
    end

    %bin tension data
    [bb,nb,sb]=bindata_line(XY,fx,bpos);
    [bc,nc,sc]=bindata_line(XY,abs(fx),bpos);
    [bv,nv,sv]=bindata2(p(:,1),v(:,1),bpos);


    subind = p(:,1)<=bpos(rl)&p(:,1)>=bpos(ll);
    cp = (XY(:,1)+XY(:,3))/2;
    subindc = cp(:,1)<=bpos(rl)&cp(:,1)>=bpos(ll);

    % store data
    stot = [stot t(ind)];  % time
    stog = [stog nanmean(diff(bv(ll:rl)')./diff(bpos(ll:rl)))]; % mean strain rate
    stof = [stof nanmean(bb(ll:rl).*nb(ll:rl))/Dy]; % mean stress (signed)
    stoa = [stoa nanmean(bc(ll:rl).*nc(ll:rl))/Dy]; % mean stress (absolute value)

    % on the timepoint we want to save draw the stress and strain profile
    if(indi==ex_indis(2))
        subplot('Position',[0.5 0.95-www*Dy/Dx_*1.18 0.4 www*Dy/Dx_*1.22])
        ax=plotyy(bpos,bb.*nb/Dy,bpos,bv);
        colorOrder = get(gca, 'ColorOrder');
        set(ax(1),'xlim',[0.5 6.5])
        set(ax(2),'xlim',[0.5 6.5])
        set(ax(1),'ylim',[0 0.005])
%         set(ax(2),'ylim',[0 0.006])
        xlabel(ax(1),'Position (\mum)') % label x-axis
        ylabel(ax(1),'Stress (nN)') % label left y-axis
        ax(1).YLabel.Color=colorOrder(1,:);
        ylabel(ax(2),'Velocity (\mum/s)') % label rigdat = [ht y-axis
    end

    % plot spatially resolved data
    if(makemovs)
        figure(h1)
        subplot(2,1,1)
        plot(p(subind,1),v(subind,1),'.')
        hold on
        plot(bpos,bv)
        plot(p(~subind,1),v(~subind,1),'.','Color',[0.75 0.75 0.75])
        hold off
        xlim([0,Dx_])
        ylim([-0.1*10^-3,1.5*10^-3])
        xlabel('x position (\mum)')
        ylabel('velocity_x (\mum/s)')
        subplot(2,1,2)
        plot(cp(subindc,1),fx(subindc),'.')
        hold on
        plot(bpos,bb)
        plot(cp(~subindc,1),fx(~subindc),'.','Color',[0.75 0.75 0.75])
        hold off
        xlim([0,Dx_])
        ylim([-0.1*10^-3,0.75*10^-3])
        ylabel('tension_x (nN)')
        xlabel('x position (\mum)')

        drawnow
        mov2(indi) = getframe(h1);
    end
    indi=indi+1;

end

%% if there was any data to store we will now display it and save it
if(makemovs)
    movie2avi(mov,[bp code '_mov_ex.avi']);
    close(h);
    movie2avi(mov2,[bp code '_data_ex.avi']);
    close(h1);
end

% only plotting the times below 80 s
sbind = stot<=80;

figure(h2);
axx=subplot('Position',[0.5 0.95-www*Dy/Dx_*2.9 0.4 0.4*Dy/Dx_*1.3]);

%% generate strain from strain rate and plot stress and strain
gam = cumtrapz(stot,stog);
[ax,p1,p2] = plotyy(stot(sbind),stof(sbind),stot(sbind),gam(sbind));

%% now edit the plot labels and colors
xlabel(ax(1),'Time (s)') % label x-axis
ylabel(ax(1),'Stress (nN)') % label left y-axis
ax(1).YLabel.Color=colorOrder(1,:);
ylabel(ax(2),'Strain') % label rigdat = [ht y-axis
% set(ax(2),'ylim',[0 0.08])
% set(ax(2),'Box','off')

%% add some stupid little lines for style points
[c, in1] = min(abs(stot-40));
axes(ax(2))
hold on
plot([4 4],[0 gam(in1)],':','Color',[0.25 0.25 0.25])
[c, in2] = min(abs(stot-10));
[c, in3] = min(abs(stot-40));
plot([10 40],[gam(in2) gam(in3)]-0.01,'--','Color',[0.25 0.25 0.25])

%% inset graph of long term tearing
pos = axx.Position;
axes('Position',[pos(1)+pos(3)*0.65 pos(2)+0.2*pos(4) pos(3)/4 pos(4)/4])
box on
plot(stot,gam,'Color',colorOrder(2,:))
xlim([0 500])
ylabel('Strain')


%% adding annotations for figure part label
annotation('textbox', [0.015 0.93 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.440 0.93 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.440 0.70 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

cd('../../../figures')
print('-depsc','-r0','figure2.eps');
cd(orig_bp)