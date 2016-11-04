function netplot_str(p,L,lf,ls,Dx,Dy,clr,clr2,str_e,str_c) 

    ncnt = ceil(L/ls)+1;
    l0 = L/(ncnt-1);

    subpL=p(mod(1:length(p),ncnt)~=0,:);
    subpR=p(mod(1:length(p),ncnt)~=1,:);
    subpL(subpL(:,1)<Dx/3&subpR(:,1)>2*Dx/3,1)=subpL(subpL(:,1)<Dx/3&subpR(:,1)>2*Dx/3,1)+Dx;
    subpL(subpL(:,2)<Dy/3&subpR(:,2)>2*Dy/3,2)=subpL(subpL(:,2)<Dy/3&subpR(:,2)>2*Dy/3,2)+Dy;
    subpR(subpR(:,1)<Dx/3&subpL(:,1)>2*Dx/3,1)=subpR(subpR(:,1)<Dx/3&subpL(:,1)>2*Dx/3,1)+Dx;
    subpR(subpR(:,2)<Dy/3&subpL(:,2)>2*Dy/3,2)=subpR(subpR(:,2)<Dy/3&subpL(:,2)>2*Dy/3,2)+Dy;
    
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
    l = l0 + l0*lf;
    str = sqrt((XY(:,3)-XY(:,1)).^2 + (XY(:,4)-XY(:,2)).^2);
    kn = (str-l)>=0;
    
%     plot(XY(:,1),XY(:,2),'.');
    
    for i =1:size(XY,1)
        if(kn(i))
            stri = max(min(length(clr),ceil(length(clr)*abs(str(i)-l)/l/str_e)),1);
            line([XY(i,1)';XY(i,3)'],[XY(i,2)';XY(i,4)'],'Color',clr(stri,:),'LineWidth',1);
        else
            stri = max(min(length(clr),ceil(length(clr)*abs(str(i)-l)/l/str_c)),1);
            line([XY(i,1)';XY(i,3)'],[XY(i,2)';XY(i,4)'],'Color',clr2(stri,:),'LineWidth',1);
        end
    end
    xlim([0 Dx]);
    ylim([0 Dy]);
end