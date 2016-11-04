function dz = activnet_pull_ode(t,z,zet,L,mu,muN,kap,xi,nu,psi,sig,Dx,Dy,Df,Dw,ncnt,lf)
% return right side of diff equation A*x'=f for node position, x
% leaves out inefficient O(~N^2) intersection step

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
    
    
    
    %% add external force at centerline and constrain edges
    if(psi>0)
        val = sig*sin(psi*t);
    elseif(psi<0)
        val = sig*round(mod(0.5+-psi*t,1)).*(round(mod(0.55+-psi*t/2,1))-0.5)*2;
    else
        val = sig;
    end
    
    subp = p(:,1)>(Df-Dw)*Dx&p(:,1)<(Df+Dw)*Dx;
    ff = 1-abs(p(subp,1)-Df*Dx)/Dw/Dx;
    if(sig<0)
        dp(subp,1)=dp(subp,1) - Dy*val.*ff/sum(ff);
    else
        dp(subp,2)=dp(subp,2) - Dy*val.*ff/sum(ff);
    end
    
    subp = p(:,1)<Dw*Dx;
    dp(subp,:)=dp(subp,:).*repmat(4*p(subp,1)/Dw/Dx-3,1,2);

    subp = p(:,1)>Dx*(1-Dw);
    dp(subp,:)=dp(subp,:).*repmat(4*abs(p(subp,1)-Dx)/Dw/Dx-3,1,2);

    subp = p(:,1)<3*Dx/4*Dw|p(:,1)>Dx-3*Dx/4*Dw;
    dp(subp,:)=0;

    %% and bring it all home
    
    dz = reshape(dp,[],1);
    
    
end