% FB Robin and Will McFadden (wmcfadden)
%
% pass in a set of points [p1x p1y;p2x p2y;p3x p3y;...etc]  and a set of
% velocities [v1x v1y;v2x v2y;v3x v3y;...etc] 
% returns barycenter p0, baryvelocity v0, velocity gradient tensor (J)
% rotation (R), strain rate tensor (E), compression rate (D), and shear
% rate tensor (S).  see http://en.wikipedia.org/wiki/Strain_rate_tensor
%
% you can also specify an optional p0 to be used as the point about which to compute the strain rates
 
function [p0,v0,J,R,E,D,S] = strainrate(p,v,p0)
    if(nargin<3)
        p0 = mean(p,1);
    end
    v0 = mean(v,1);
    A = [];
    dv = [];
    for i=1:size(p,1)
        A = [A; p(i,1)-p0(1,1) p(i,2)-p0(1,2) 0 0; 0 0 p(i,1)-p0(1,1) p(i,2)-p0(1,2)];
        dv = [dv; v(i,1)-v0(1,1); v(i,2)-v0(1,2)];
    end
    J = (A'*A)\A'*dv;
    J = reshape(J,2,2)';
    R = (J-J')/2;
    E = (J+J')/2;
    D = trace(E)/2;
    S = E-D*eye(2);
end