function netplot(p,L,lf,ls,Dx,Dy,clr)    
    ncnt = ceil(L/ls)+1;
    l0 = L/(ncnt-1);

    subpL = p;
    subpR = p;
    subpL=subpL(mod(1:length(subpL),ncnt)~=0,:);
    subpR=subpR(mod(1:length(subpR),ncnt)~=1,:);
    subpL(subpL(:,1)<Dx/4&subpR(:,1)>3*Dx/4,1)=subpL(subpL(:,1)<Dx/4&subpR(:,1)>3*Dx/4,1)+Dx;
    subpL(subpL(:,2)<Dy/4&subpR(:,2)>2*Dy/4,2)=subpL(subpL(:,2)<Dy/4&subpR(:,2)>2*Dy/4,2)+Dy;
    subpR(subpR(:,1)<Dx/4&subpL(:,1)>3*Dx/4,1)=subpR(subpR(:,1)<Dx/4&subpL(:,1)>3*Dx/4,1)+Dx;
    subpR(subpR(:,2)<Dy/4&subpL(:,2)>2*Dy/4,2)=subpR(subpR(:,2)<Dy/4&subpL(:,2)>2*Dy/4,2)+Dy;
    
    subv = subpR-subpL;
    subv = subv./repmat(sqrt(subv(:,1).^2+subv(:,2).^2),1,2);
    subpL = subpL - l0*lf/2*subv;
    subpR = subpR + l0*lf/2*subv;
    
    
    XY = [subpL subpR];
    subXY = XY(:,1)>Dx|XY(:,2)>Dy|XY(:,3)>Dx|XY(:,4)>Dy;
    
    extXY = XY(subXY, :);
    
    tsub = extXY(:,1)>Dx;
    extXY(tsub,:)=extXY(tsub,:)-repmat([Dx 0 Dx 0],sum(tsub),1);
    tsub = extXY(:,3)>Dx;
    extXY(tsub,:)=extXY(tsub,:)-repmat([Dx 0 Dx 0],sum(tsub),1);
    tsub = extXY(:,2)>Dy;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 Dy 0 Dy],sum(tsub),1);
    tsub = extXY(:,4)>Dy;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 Dy 0 Dy],sum(tsub),1);
    
    XY = [XY; extXY];
    line([XY(:,1)';XY(:,3)'],[XY(:,2)';XY(:,4)'],'Color',clr,'LineWidth',2);
    xlim([0 Dx]);
    ylim([0 Dy]);
end