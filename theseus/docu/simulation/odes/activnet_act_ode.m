function dz = activnet_act_ode(t,z,zet,L,mu,muN,kap,xi,nu,psi,sig,Dx,Dy,Df,Dw,ncnt,lf)
% return right side of diff equation A*x'=f for node position, x


    %% compute intrafilament forces    
    l0 = L/(ncnt-1);
    p = reshape(z,[],2);
    p = [mod(p(:,1),Dx),mod(p(:,2),Dy)];
    dp = zeros(size(p));
    for n=1:ncnt:length(p)
        va_orth=[0 0];
        va = [0 0];
        la = 0;
        for i=0:ncnt-2
            vb = mydiff(p(n+i,:),p(n+i+1,:),Dx,Dy);
            lb = sqrt(vb*vb');
            gam = (lb-l0)/l0;
            f = mu*vb/lb*gam;
            if(mu<0)
                f = -f*(1+(muN-1)*double(gam>0));
            end
            dp(n+i,:) = dp(n+i,:) + f;
            dp(n+i+1,:) = dp(n+i+1,:) - f;
            vb_orth = [-vb(2) vb(1)];
            if(i>0)
                if(va_orth*vb'>0)
                    va_orth = -va_orth;
                end
                if(vb_orth*va'<0)
                    vb_orth = -vb_orth;
                end
                tor = kap/l0^2*acos(max(min(va*vb'/la/lb,1),0));
                dp(n+i-1,:)=dp(n+i-1,:)+tor*va_orth/la;
                dp(n+i,:)=dp(n+i,:)-tor*va_orth/la;
                dp(n+i+1,:)=dp(n+i+1,:)+tor*vb_orth/lb;
                dp(n+i,:)=dp(n+i,:)-tor*vb_orth/lb;
            end
            va = vb;
            va_orth = vb_orth;
            la = lb;
        end
    end
    
    
    %% add active force from motors at crosslinking points
    if(~isempty(nu))
        indL = 1:length(p);
        indL = indL(mod(indL,ncnt)~=0);

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

        indL = [indL indL(subXY)];

        g = lineSegmentGrid(indL,XY,Dx,Dy,l0);

        f = min(1,max(0,(g-lf/2)/(1-lf)));
        for ind=1:size(g,1)
            i = g(ind,3);
            j = g(ind,4);
            
            vm = mydiff(p(j,:),p(j+1,:),Dx,Dy);
            lm = sqrt(vm*vm');
            
            edg = 1;
            
            if(g(ind,1)<lf)
                edg = edg*g(ind,1)/lf;
            elseif((1-g(ind,1))<lf)
                edg = edg*(1-g(ind,1))/lf;
            end
            
            if(g(ind,2)<lf)
                edg = edg*g(ind,2)/lf;
            elseif((1-g(ind,2))<lf)
                edg = edg*(1-g(ind,2))/lf;
            end
            mul = 1;
            if(psi>0)
                mul = double(g(ind,5)>=psi*abs(-Dx*Df));
            end
            tnu = nu(ceil(i/ncnt),ceil(j/ncnt))*mul;
            dp(i:i+1,:) = dp(i:i+1,:) + edg*tnu/lm*[vm*(1-f(ind,1));vm*f(ind,1)];
            dp(j:j+1,:) = dp(j:j+1,:) - edg*tnu/lm*[vm*(1-f(ind,2));vm*f(ind,2)];
            
        end
    end
    
    
    %% and bring it all home
    
    dz = reshape(dp,[],1);
    
    
end