bp = '~/Documents/MATLAB/activnet/data/examples/fig6/';
code = 'vimlpcks';%gcqbbcyr
cd(bp)

makemovs = 0;

%% load param file and decipher params
fid = fopen([bp code '_scr.txt']);
C = textscan(fid, '%s','delimiter', '\n');
fclose(fid);
pare = strsplit(C{1}{9}, '>');
paree = strsplit(pare{1}, ' ');
paree = {paree{2:end}};
zet=str2num(paree{2});L=str2num(paree{3});mu=str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6}); 
xi=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
r=str2num(paree{11});sig=str2num(paree{12});Dx=str2num(paree{13});Dy=str2num(paree{14});Df=str2num(paree{15});
Dw=str2num(paree{16});ls=str2num(paree{17});lf=str2num(paree{18});


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
stot = [];
stog = [];
stoe = [];
stoc = [];
gp = [];
sp = [];

%% setup timepoints and space points to measure
inds = 1:10:size(zt,1);
inds = inds(2:end);
ex_indis = [1 50 400];
bpos = linspace(0,Dx,51);
bpos = bpos(1:end-1)+bpos(2)/2;
ll = 10;
rl = 21;

%% for loop over timepoints to display
h2=figure;
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
Dx_ = Dx;
exw = 0.3;
cdmn = 0.5;
edmn = 0.05;
lpr = 10;
rpr = 90;

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
    [XY,sx,sy,str]=get_str(p,L,lf,ls,Dx,Dy);   

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
    cpx = (XY(:,1)+XY(:,3))/2;
    subindc = cpx(:,1)<=bpos(rl)&cpx(:,1)>=bpos(ll);
    gp = [gp; bv'];
    sp = [sp; bb];
    % store data
    stot = [stot t(ind)];
    stog = [stog nanmean(diff(bv(ll:rl)')./diff(bpos(ll:rl)))];
    stoe = [stoe mean(str(str>0))];
    stoc = [stoc mean(str(str<0))];
   
%     if(indi==ex_indis(2))
%         subplot('Position',[0.05 0.925-exw*Dy/Dx_-0.2 0.2 0.2])
%         ax=plotyy(bpos,bb.*nb/Dy,bpos(1:length(bv)),bv*10);
%         colorOrder = get(gca, 'ColorOrder');
%         set(ax(1),'xlim',[0 Dx])
%         set(ax(2),'xlim',[0 Dx])
%         xlabel(ax(1),'Position (\mum)') % label x-axis
%         ylabel(ax(1),'Stress (nN)') % label left y-axis
%         ax(1).YLabel.Color=colorOrder(1,:);
%         ylabel(ax(2),'Velocity (\mum/s)') % label rigdat = [ht y-axis
%     end
    
    % plot spatially resolved data 
    if(makemovs)
        figure(h1)
        subplot(2,1,1)
        plot(p(~subind,1),v(~subind,1),'.')
        hold on
        plot(p(subind,1),v(subind,1),'.')
        plot(bpos,bv)
        hold off
        xlim([0,Dx_])
        ylim([-0.01,0.01])
        xlabel('x position (\mum)')
        ylabel('velocity_x (\mum/s)')
        subplot(2,1,2)
        plot(cpx(~subindc,1),fx(~subindc),'.')
        hold on
        plot(cpx(subindc,1),fx(subindc),'.')
        plot(bpos,bb)
        hold off
        xlim([0,Dx_])
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
figure(h2)
subplot('Position',[0.125 0.95-exw*Dy/Dx_ exw exw*Dy/Dx_])
[ax,p1,p2]=plotyy(bpos,mean(gp),bpos(5:end),mean(sp(:,5:end)));
ylim([-0.003,0.003])
xlabel('position (\mum)')
ylabel('velocity (\mum/s)')
set(ax(1),'YTick',[-0.003,-0.0015,0,0.0015,0.003],'YTickLabel',[-0.003,-0.0015,0,0.0015,0.003])
axes(ax(2))
ylim([0,0.005])
set(ax(2),'YTick',[0,0.001,0.002,0.003,0.004,0.005],'YTickLabel',[0,0.001,0.002,0.003,0.004,0.005])
ylabel('stress (nN/\mum)')

title('\tau_r=10000 s')


