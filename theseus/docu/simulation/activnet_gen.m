function p = activnet_gen(zet,L,mu,kap,lc,xi,ups,phi,psi,r,sig,Dx,Dy,Df,Dw,ls,lf,tinc,tfin,nonlin)
% generates an active network simulation and prints node positions
% at time steps.  Parameters are defined as follows:
%
%   zet - medium viscosity
%   L - length of the filament
%   mu - compressional modulus of the filament
%   kap - bending modulus of a filament if ls<L
%   lc - average distance between filament overlaps
%   xi - frictional resistance between two overlapping segments
%   ups - motor force at filament overlaps
%   phi - fraction of overlaps that receive a motor force
%   psi - spatial variation in motor force (see below)
%   sig - applied force
%   Dx - x-dimension of domain
%   Dy - y-dimension of domain
%   Df - 
%   Dw - width of window in x-dimension where forces/constraints are applied
%   ls - length of filament segments
%   lf - length of force falloff at end of filament (for continuous forces)
%   tinc - time increment to return solutions
%   tfin - end time of simulation
%   nonlin - nonlinear factor by which to make filament stiffer by extension
%
%   See also SUM, PLUS.    


%% this mess just ensures that any string input is converted to numbers
    if(ischar(zet)); zet = str2num(zet); end;
    if(ischar(L)); L = str2num(L); end;
    if(ischar(mu)); mu = str2num(mu); end;
    if(ischar(kap)); kap = str2num(kap); end;
    if(ischar(lc)); lc = str2num(lc); end;
    if(ischar(xi)); xi = str2num(xi); end;
    if(ischar(ups)); ups = str2num(ups); end;
    if(ischar(phi)); phi = str2num(phi); end;
    if(ischar(psi)); psi = str2num(psi); end;
    if(ischar(r)); r = str2num(r); end ;
    if(ischar(sig)); sig = str2num(sig); end;
    if(ischar(Dx)); Dx = str2num(Dx); end;
    if(ischar(Dy)); Dy = str2num(Dy); end;
    if(ischar(Df)); Df = str2num(Df); end;
    if(ischar(Dw)); Dw = str2num(Dw); end;
    if(ischar(ls)); ls = str2num(ls); end;
    if(ischar(lf)); lf = str2num(lf); end;
    if(ischar(tinc)); tinc = str2num(tinc); end;
    if(ischar(tfin)); tfin = str2num(tfin); end;
    if(ischar(nonlin)); nonlin = str2num(nonlin); end;
    Dp = 1;
    if(Df<0);Df=abs(Df);Dp = Df; end;
%     rng(abs(nonlin));
    
    %% use inputs to calculate number of filaments to add
    ncnt = ceil(L/ls)+1;
    N = floor(2*Dx*Dy/lc/L);
    
    nu=[];
    if(ups>0)
        nu = ups*double(rand(N,N)<phi).*(ones(N,N)-eye(N,N));
        nu = (nu+nu')/2;
    end
    
    %% initialize network
    
    p = zeros(N*ncnt,2);
    for i=1:N
        p((i-1)*ncnt+1,:) = [Dp*Dx*rand Dy*rand];
        thet = rand*2*pi;
        for j = 2:ncnt
            p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1.0)*[cos(thet) sin(thet)];
        end
    end
    if(nonlin<0)
        p = zeros(N*ncnt,2);
        for i=1:N
            p((i-1)*ncnt+1,:) = [Dx*(0.2+0.6*rand) Dy*(0.2+0.6*rand)];
            thet = rand*2*pi;
            for j = 2:ncnt
                p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1.0)*[cos(thet) sin(thet)];
            end
        end
        nonlin=-nonlin;
    end
    
    muN = nonlin;
    
    p = [mod(p(:,1),Dx),mod(p(:,2),Dy)];
    
    fileID = 1;
    %% solve ode
    z0 = reshape(p,1,[]);

    if(tinc>0.05/2/r)
        tinc = 0.05/2/r;
    end
    tt = 0:tinc:tfin;
    
    fprintf(fileID,'%.3f',0);
    for i=1:length(z0)
        fprintf(fileID,' %.4f',z0(i));
    end
    fprintf(fileID,'\n');
    
    activnet(N,tt,z0,zet,L,mu,muN,kap,xi,nu,psi,sig,Dx,Dy,Df,Dw,Dp,ncnt,lf,r,tinc,fileID);
    

end
