    
%% this mess just ensures that any string input is converted to numbers
zet = 0.005;
L = 4;
mu = 0.1;
kap = 0;
xi = 1;
ups = 1;
phi = 1;
psi = 1;
r = 0;
sig = 0;
Dx = 10;
Dy = 10;
Df = 0.5;
Dw = 0.1;
ls = 4;
lf = 0.05;
tinc = 0.01;
tfin = 10;
seed = 100;

rng(seed);

%% use inputs to calculate number of filaments to add
ncnt = ceil(L/ls)+1;
N = 2;

nu=[];
if(ups>0)
    nu = ups*double(rand(N,N)<phi).*(ones(N,N)-eye(N,N));
    nu = (nu+nu')/2;
end

%% initialize network

p = zeros(N*ncnt,2);
i=1;
p((i-1)*ncnt+1,:) = [2.5 3];
thet = 0;
for j = 2:ncnt
    p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1.0)*[cos(thet) sin(thet)];
end
i=2;
p((i-1)*ncnt+1,:) = [3 2.5];
thet = pi/2;
for j = 2:ncnt
    p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1.0)*[cos(thet) sin(thet)];
end

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

activnet(N,tt,z0,zet,L,mu,kap,xi,nu,psi,sig,Dx,Dy,Df,Dw,ncnt,lf,r,tinc,fileID);
