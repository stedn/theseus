function v = mydiff(p1,p2,Dx,Dy)
% specialized vector differences in the periodic domain
    v = [mysub(p1(:,1),p2(:,1),Dx) mysub(p1(:,2),p2(:,2),Dy)];
end