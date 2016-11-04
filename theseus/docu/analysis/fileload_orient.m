% bp = '/Users/wmcfadden/xlrelax_all/';
% code = 'vghjadsd';
cd(bp)
fid = fopen([bp code '_scr.txt']);
C = textscan(fid, '%s','delimiter', '\n');
pare = strsplit(C{1}{9}, '>');
paree = strsplit(pare{1}, ' ');
paree = {paree{2:end}};
zet=str2num(paree{2});L=str2num(paree{3});mu=str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6}); 
del=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
r=str2num(paree{11});sig=str2num(paree{12});D=str2num(paree{13});Df=str2num(paree{14});ls=str2num(paree{15});lf=str2num(paree{16});
% zet=str2num(paree{2});L=str2num(paree{3});mu=str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6}); 
% del=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
% %         r=str2num(paree{10});sig=str2num(paree{11});D=str2num(paree{12});Df=str2num(paree{13});ls=str2num(paree{14});lf=str2num(paree{15});
% r=str2num(paree{11});sig=str2num(paree{12});D=str2num(paree{13});Df=str2num(paree{14});ls=str2num(paree{15});lf=str2num(paree{16});
%         sig=str2num(paree{10});D=str2num(paree{11});Df=str2num(paree{12});ls=str2num(paree{13});lf=str2num(paree{14});
fclose(fid);
if(L/lc>6)
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
    xt = zt(:,1:end/2);
    yt = zt(:,end/2+1:end);
    %             dy = yt-repmat(yt(1,:),size(yt,1),1);
    if(length(t)>1)
        coff = abs(2*median(yt(2,:)-yt(1,:)));
        dy = 0*yt;
        for k=2:size(yt,1)
            dyc = yt(k,:)-yt(k-1,:);
            subs = abs(dyc)<coff;
            dy(k,subs)=dy(k-1,subs)+dyc(subs);
            dy(k,~subs)=dy(k-1,~subs);
        end

        nbins = 30;
        brng = linspace(0,2*D,nbins+1)';
        sp = brng(1:end-1)+brng(2)/2;
        sp(sp>D)=2*D-sp(sp>D);

        ss = [];
        tt = [];
        for k=1:size(xt,1)
            if(sum(isnan(xt(k,:)))==0)
                sv = bindata(dy(k,:),xt(k,:),brng);
                subs = sp>Df*D&sp<D*(1-Df);
                tt=[tt; t(k)];
                ss=[ss; -nanmean(sv(subs)./sp(subs))];
            end
        end

%         clear mov
%         figure('Position', [50, 100, 600, 300]);
        lst = size(zt,1);
        trp = repmat((1:lst)'/lst,1,3);
        cc = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*copper(lst);
        temp=hot(2*lst);
        cc2 = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*temp(1:lst,:);
        edges = {linspace(0.65,1.35,50),linspace(-pi/2,pi/2,20)}  ;      
        indi = 1;
        stohist = [];
        strs = [];
        stostr = [];
        stopos = [];
        stoang = [];
        stot = [];
        stostd = [];
        for ind = 1:size(zt,1)
            p = reshape(zt(ind,:),[],2);
            p = [mod(p(:,1),2*D),mod(p(:,2),D)];
        %     whitebg('black')
        %     set(gcf,'Color',[0 0 0])
        %     set(gcf,'InvertHardcopy','off')

        %     netplot_str(p,L,lf,ls,D,cc,cc2,0.25);
            ncnt = ceil(L/ls)+1;
            l0 = L/(ncnt-1);
            subpL = p;
            subpR = p;
            subpL=subpL(mod(1:length(subpL),ncnt)~=0,:);
            subpR=subpR(mod(1:length(subpR),ncnt)~=1,:);
            subpL(subpL(:,1)<D/2&subpR(:,1)>3*D/2,1)=subpL(subpL(:,1)<D/2&subpR(:,1)>3*D/2,1)+2*D;
            subpL(subpL(:,2)<D/3&subpR(:,2)>2*D/3,2)=subpL(subpL(:,2)<D/3&subpR(:,2)>2*D/3,2)+D;
            subpR(subpR(:,1)<D/2&subpL(:,1)>3*D/2,1)=subpR(subpR(:,1)<D/2&subpL(:,1)>3*D/2,1)+2*D;
            subpR(subpR(:,2)<D/3&subpL(:,2)>2*D/3,2)=subpR(subpR(:,2)<D/3&subpL(:,2)>2*D/3,2)+D;
            subv = subpR-subpL;
            subv = subv./repmat(sqrt(subv(:,1).^2+subv(:,2).^2),1,2);
            subpL = subpL;
            subpR = subpR;
            XY = [subpL subpR];
            str = sqrt((XY(:,3)-XY(:,1)).^2 + (XY(:,4)-XY(:,2)).^2);
            angs = atan((XY(:,4)-XY(:,2))./(XY(:,3)-XY(:,1)));
            subind = XY(:,1)<D*(1-Df)&XY(:,3)<D*(1-Df)&XY(:,1)>D*Df&XY(:,3)>D*Df;
            pos = (XY(:,1)+XY(:,3))/2;

%             strs = [strs; mean(abs(str(subind)-L)/L)];
            stopos = [stopos; pos(subind)];
            stostr = [stostr; mean((str(subind)-L).^2)];
            stoang = [stoang; circ_mean(2*angs(subind))/2];
            stot = [stot; tt(ind)*ones(size(angs(subind)))];
            x = angs(subind);
            y = (str(subind)-L)/L;
            n=20;
            bin = linspace(min(x), max(x), n+1);
            ind = sum(bsxfun(@minus, x, ones(length(x),1)*bin)>=0,2);

            m = NaN(n,1);
            e = NaN(n,1);
            s = NaN(n,1);
            for i = 1:n
                m(i) = mean(y(ind==i));   % Mean value over the bin
                e(i) = std(y(ind==i));    % Standard deviation
                s(i) = sum((ind==i));    % sum
            end
            strs = [strs; sqrt(2*mean(y.^2))];
            stostd = [stostd; sum(e.*s)./sum(s)];
%             stostd = [stostd; circ_std(2*angs(subind))/2];  % mean of standard deviations
        %     plot(x,y,'.');
        %     hold on
        %     u = (bin(1:end-1)+bin(2:end))/2;
        %     errorbar(u,m,e,'k');

%             h = histogram(angs(subind),edges{2},'Normalization','probability');
%             stohist = [stohist; h.Values];
%             h = histogram(str(subind)/L,edges{1},'Normalization','probability');
%             stohist = [stohist; h.Values];
        %     hist3([str(subind)/L,angs(subind)],'Edges',edges);
        %     [N,C] = hist3([str(subind)/L,angs(subind)],'Edges',edges);
        %     imagesc(C{2},C{1},N)
%             scatter(angs(subind),str(subind)/L, 5, str(subind)/L);
%             colormap([flipud(cc2);cc])
%             set(gca,'CLim',[0.85 1.95]);
%             ylim([0.85,1.15]);
%             mov(indi) = getframe;
%             clf
%             indi = indi +1;
        end
%         movie2avi(mov,[code '_dist_mov.avi']);

        % nbins = 30;
        % rs = linspace(0,1,nbins+1);
        % rs = rs(1:end-1);
        % sp = D*(rs+rs(2)/2);
        % 
        % figure;
        % xt = zt(:,1:end/2);
        % yt = zt(:,end/2+1:end);
        % for k=1:size(xt,1)
        %     sv = [];
        %     for r=rs
        %         subx = xt(:,max(xt(:,1:end/2))<(r+0.05)*D&min(xt(:,1:end/2))>r*D);
        %         suby = yt(:,max(xt(:,1:end/2))<(r+0.05)*D&min(xt(:,1:end/2))>r*D);
        %         sv = [sv nanmean(suby(k,:)-(suby(1,:)))];
        %     end
        %     ss(k)=-mean(sv./sp);
        %     sr(k)=std(sp*ss(k)-sv);
        % end
        % xsi0 = 1;
        % gam0 = 
        % fito = fit(t,ss','(1-exp(-xsi*x))*gam','StartPoint', [xsi0, gam0]);
        n=10;
        subss= ss(1:5:size(zt,1));
        substr= strs(1:5:size(zt,1));
        substd= stostd(1:5:size(zt,1));
        subtt=tt(1:5:size(zt,1));
        kr = lc.^2./zet./del./(L-2*lc).^2*4*pi;
        x = subtt;
        y = ((stostd));
        bin = linspace(min(x), max(x), n+1);
        ind = sum(bsxfun(@minus, x, ones(length(x),1)*bin)>=0,2);
        m = NaN(n,1);
        e = NaN(n,1);
        s = NaN(n,1);
        for i = 1:n
            m(i) = mean(y(ind==i));   % Mean value over the bin
            e(i) = std(y(ind==i));    % Standard deviation
            s(i) = sum(y(ind==i));    % Standard deviation
        end
        mm = m;
        y = ((substr));
        bin = linspace(min(x), max(x), n+1);
        ind = sum(bsxfun(@minus, x, ones(length(x),1)*bin)>=0,2);
        m = NaN(n,1);
        e = NaN(n,1);
        s = NaN(n,1);
        for i = 1:n
            m(i) = mean(y(ind==i));   % Mean value over the bin
            e(i) = std(y(ind==i));    % Standard deviation
            s(i) = sum(y(ind==i));    % Standard deviation
        end
        mmm = m;
        y = ((subtt));
        bin = linspace(min(x), max(x), n+1);
        ind = sum(bsxfun(@minus, x, ones(length(x),1)*bin)>=0,2);
        m = NaN(n,1);
        e = NaN(n,1);
        s = NaN(n,1);
        for i = 1:n
            m(i) = mean(y(ind==i));   % Mean value over the bin
            e(i) = std(y(ind==i));    % Standard deviation
            s(i) = sum(y(ind==i));    % Standard deviation
        end
        tt = m;
        y = ((stostr));
        bin = linspace(min(x), max(x), n+1);
        ind = sum(bsxfun(@minus, x, ones(length(x),1)*bin)>=0,2);
        m = NaN(n,1);
        e = NaN(n,1);
        s = NaN(n,1);
        for i = 1:n
            m(i) = mean(y(ind==i));   % Mean value over the bin
            e(i) = std(y(ind==i));    % Standard deviation
            s(i) = sum(y(ind==i));    % Standard deviation
        end
        mmmm = m;
        y = ((subss));
        bin = linspace(min(x), max(x), n+1);
        ind = sum(bsxfun(@minus, x, ones(length(x),1)*bin)>=0,2);
        m = NaN(n,1);
        e = NaN(n,1);
        s = NaN(n,1);
        for i = 1:n
            m(i) = mean(y(ind==i));   % Mean value over the bin
            e(i) = std(y(ind==i));    % Standard deviation
            s(i) = sum(y(ind==i));    % Standard deviation
        end
        
%         mmmm is the total energy
%         mmm is the stretch
%         mm is the stretch dispersion
%         m is the strain
        alltt = [alltt tt];
        allenrg = [allenrg mm];
        alldisp = [alldisp diff(mm./max(mmm))./diff(tt)];
        allstrain = [allstrain (m/mmm(end))];
        allstrun = [allstrun m];
        allvisc = [allvisc mmm];
        alltest = [alltest diff(log(m))./diff(log(tt))];
        allangs = [allangs; stoang];
        details = [details [zet;L;mu;kap;lc;del;sig]];
%         subplot(2,1,1)
%         loglog(subtt,subss/substr(end),'DisplayName',['\mu = ' num2str(mu) ', \xi = ' num2str(del*zet)])
%         hold on
%         loglog(subtt,substr/substr(end),'--')
%         subplot(2,1,2)
%         semilogx(subtt,substd./substr,'DisplayName',['\mu = ' num2str(mu) ', \xi = ' num2str(del*zet)])
%         hold on
%         drawnow
    end
end