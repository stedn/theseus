function [XY, sx, sy, str] = get_str(p,L,lf,ls,Dx,Dy)    
    ncnt = ceil(L/ls)+1;
    l0 = L/(ncnt-1);

    subpL=p(mod(1:length(p),ncnt)==1,:);
    subpR=p(mod(1:length(p),ncnt)==0,:);
    subpL(subpL(:,1)<Dx/3&subpR(:,1)>2*Dx/3,1)=subpL(subpL(:,1)<Dx/3&subpR(:,1)>2*Dx/3,1)+Dx;
    subpL(subpL(:,2)<Dy/3&subpR(:,2)>2*Dy/3,2)=subpL(subpL(:,2)<Dy/3&subpR(:,2)>2*Dy/3,2)+Dy;
    subpR(subpR(:,1)<Dx/3&subpL(:,1)>2*Dx/3,1)=subpR(subpR(:,1)<Dx/3&subpL(:,1)>2*Dx/3,1)+Dx;
    subpR(subpR(:,2)<Dy/3&subpL(:,2)>2*Dy/3,2)=subpR(subpR(:,2)<Dy/3&subpL(:,2)>2*Dy/3,2)+Dy;
    
    XY = [subpL subpR];
    
    ang = atan2((XY(:,4)-XY(:,2)),(XY(:,3)-XY(:,1)));
    str = (sqrt((XY(:,3)-XY(:,1)).^2 + (XY(:,4)-XY(:,2)).^2)-l0)/l0;
    sx = str.*abs(cos(ang));
    sy = str.*abs(sin(ang));
    
end