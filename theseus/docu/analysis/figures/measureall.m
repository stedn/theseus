cd(bp)

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
Dw=str2num(paree{16});ls=str2num(paree{17});lf=str2num(paree{18});nonl=str2num(paree{21});


%% load simulation data
t = [];
A = [];
if(mu<0)
    A = importdata([bp code '_out.txt']);
    if(isstruct(A))
        A = A.data;
        if(size(A,1)==1)
            imp2 = importdata([bp code '_out.txt'],' ',9);
            if(isfield(imp2,'data'))
                A = [A;imp2.data];
            end
        end
        t = A(:,1);
        zt = A(:,2:end);
   end
end
if(isempty(t))
    t=0;
    zt = [0 0 1 1];
end
%% store initial positions and initialize measurements
op = reshape(zt(1,:),[],2);
tl=0;
stof = [];
stog = [];
stoa = [];
stot = [];
stoe = [];
stoc = [];
stofe = [];
stofc = [];
stow = [];

   
%% setup timepoints and space points to measure
inds = 1:10:size(zt,1);
inds = inds(2:end);
bpos = linspace(0,Dx,bns);
bpos = bpos(1:end-1)+bpos(2)/2;
lpr = 10;
rpr = 90;

%% for loop over timepoints to measure


for ind = inds
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),Dx),mod(p(:,2),Dy)];
    
    dp = (p-op);
    op = p;
    
    % remove data if has moved farther than realistically possible 
    % these events are due to crossing domain boundary or recycling
    jumpcut = max(10*median(abs(dp(1:1:end,1))),0.0005);
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
    cpx = (XY(:,1)+XY(:,3))/2;
    
    % compute filament tension
    fstr = mu*str;
    fx = mu*sx;
    fy = mu*sy;
    if(mu<0)
        fx = -fx.*(1+99*double(sx>0));
        fy = -fy.*(1+99*double(sy>0));
        fstr = -fstr.*(1+99*double(str>0));
    end

    %bin tension data
    [bb,nb,sb]=bindata_line(XY,fx,bpos);
    [bc,nc,sc]=bindata_line(XY,abs(fx),bpos);
    [bv,nv,sv]=bindata2(p(:,1),v(:,1),bpos);
    
    subind = p(:,1)<=bpos(rl)&p(:,1)>=bpos(ll);
    % store data
    stot = [stot t(ind)];
    stog = [stog nanmean(diff(bv(ll:rl)')./diff(bpos(ll:rl)))];
    stow = [stow (prctile(cpx,rpr)-prctile(cpx,lpr))];
    stof = [stof nanmean(bb(ll:rl).*nb(ll:rl))/Dy];
    stoa = [stoa nanmean(bc(ll:rl).*nc(ll:rl))/Dy];
    stoe = [stoe mean(str(str>0))];
    stoc = [stoc mean(str(str<0))];
    stofe = [stofe sum(fstr(str>0))];
    stofc = [stofc sum(fstr(str<0))];
end


%% if there was any data to store we will now display it and save it
if(length(stof)>2)
    
    % saving data, must extend sto data to fit with all data
    if(~isempty(allg))
        stog = [stog zeros(1,size(allg,2)-length(stog))];
        stof = [stof zeros(1,size(allg,2)-length(stof))];
        stot = [stot zeros(1,size(allg,2)-length(stot))];
        stoc = [stoc zeros(1,size(allc,2)-length(stoc))];
        stoe = [stoe zeros(1,size(alle,2)-length(stoe))];
%        stofc = [stofc zeros(1,size(allfc,2)-length(stofc))];
%        stofe = [stofe zeros(1,size(allfe,2)-length(stofe))];
        stoa = [stoa zeros(1,size(alla,2)-length(stoa))];
        stow = [stow zeros(1,size(allw,2)-length(stow))];
        allg = [allg zeros(size(allg,1),length(stog)-size(allg,2))];
        allf = [allf zeros(size(allf,1),length(stof)-size(allf,2))];
        allt = [allt zeros(size(allt,1),length(stot)-size(allt,2))];
        allc = [allc zeros(size(allc,1),length(stoc)-size(allc,2))];
        alle = [alle zeros(size(alle,1),length(stoe)-size(alle,2))];
%        allfc = [allfc zeros(size(allfc,1),length(stofc)-size(allfc,2))];
%        allfe = [allfe zeros(size(allfe,1),length(stofe)-size(allfe,2))];
        alla = [alla zeros(size(alla,1),length(stoa)-size(alla,2))];
        allw = [allw zeros(size(allw,1),length(stow)-size(allw,2))];
    end
    % concatenate data with all matrices
    allg = [allg; stog];
    allf = [allf; stof];
    allt = [allt; stot];
    allc = [allc; stoc];
    alle = [alle; stoe];
    alla = [alla; stoa];
    allw = [allw; stow];
    allp = [allp; zet L mu kap lc xi ups phi psi r sig Dx Dy Df Dw nonl];
    alln = {alln{:} code};
end
