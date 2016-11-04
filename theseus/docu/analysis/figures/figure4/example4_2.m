code = 'ekwfpndo'

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
t = A(:,1);
zt = A(:,2:end);

%% store initial positions and initial measurements (all 0)
op = reshape(zt(1,:),[],2);
tl=0;
stof2 = 0;
stog2 = 0;
stoa2 = 0;
stot2 = 0;


%% setup timepoints and space points to measure
[c, lind] = min(abs(t-380));
inds = 1:20:lind;
inds = inds(2:end);
ex_indis = [1 19];
bpos = linspace(0,Dx,51);
bpos = bpos(1:end-1)+bpos(2)/2;
ll = 4;
rl = 16;

if(makemovs)
    clear mov mov2
    h1 = figure;
    h = figure('Position', [50, 100, 100+600*Dx/Dy, 600]);
end
lst = size(zt,1);
trp = repmat((1:lst)'/lst,1,3);
temp=flipud(winter(lst));
temp2 = bone(2*lst);
cc = (1-trp.^2).*temp2(lst+1:end,:)+(trp.^2).*temp(1:lst,:);
temp=hot(2*lst);
cc2 = (1-trp.^2).*temp2(lst+1:end,:)+(trp.^2).*temp(1:lst,:);
indi = 1;

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
        colorbar('westoutside')
        drawnow
        mov(indi) = getframe(h);
    end
    if(indi==ex_indis(1))
        subplot('Position',[0.05 0.915-0.2*Dy/Dx_*2 0.2 0.2*Dy/Dx_])
        pat=patch([Dx_-Dx*Dw Dx_-Dx*Dw Dx_ Dx_],[0 Dy Dy 0],[.7 .5 0]);
        set(pat,'FaceAlpha',0.25,'EdgeColor','none');
        netplot_str(p,L,lf,ls,Dx,Dy,cc,cc2,edmn,cdmn);
        ylabel(['\tau_r = ' num2str(1/r) ' s'])
        xlim([0 Dx_])
        set(gca,'xtick',[],'ytick',[],'box','on')
        xlabel('\gamma = 0')

    end
    if(indi==ex_indis(2))
        subplot('Position',[0.275 0.915-0.2*Dy/Dx_*2 0.2 0.2*Dy/Dx_])
        pat=patch([Dx_-Dx*Dw Dx_-Dx*Dw Dx_ Dx_],[0 Dy Dy 0],[.7 .5 0]);
        set(pat,'FaceAlpha',0.25,'EdgeColor','none');
        netplot_str(p,L,lf,ls,Dx,Dy,cc,cc2,edmn,cdmn);
        xlim([0 Dx_])
        set(gca,'xtick',[],'ytick',[],'box','on')
        xlabel('\gamma = 0.4')

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
    bv=bindata(v(:,1),p(:,1),bpos);


    subind = p(:,1)<=bpos(rl)&p(:,1)>=bpos(ll);
    cp = (XY(:,1)+XY(:,3))/2;
    subindc = cp(:,1)<=bpos(rl)&cp(:,1)>=bpos(ll);

    % store data
    stot2 = [stot2 t(ind)];
    stog2 = [stog2 nanmean(v(subind,1)./(p(subind,1)-Dx*Dw))];
    stof2 = [stof2 nanmean(bb(ll:rl).*nb(ll:rl))/Dy];
    stoa2 = [stoa2 nanmean(bc(ll:rl).*nc(ll:rl))/Dy];



    % plot spatially resolved data
    if(makemovs)
        figure(h1)
        subplot(2,1,1)
        plot(p(~subind,1),v(~subind,1),'.')
        hold on
        plot(p(subind,1),v(subind,1),'.')
        plot(bpos(1:end-1),bv)
        hold off
        xlim([0,Dx_])
        ylim([-0.1*10^-3,0.5*10^-2])
        xlabel('x position (\mum)')
        ylabel('velocity_x (\mum/s)')

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
