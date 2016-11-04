p = rand([400,2]);
v = [p(:,1) p(:,2)];
[p0,v0,J,R,E,D,S] = strainrate(p,v);

figure;
pc = [p(:,1)-p0(1) p(:,2)-p0(2)];
quiver(pc(:,1),pc(:,2),v(:,1),v(:,2));


figure;
vc = [v(:,1)-v0(1) v(:,2)-v0(2)];
quiver(pc(:,1),pc(:,2),vc(:,1),vc(:,2));
hold on;
quiver(0,0,v0(1),v0(2),'r')

figure;
vj = [];
for i=1:length(vc)
    vj = [vj J*pc(i,:)'];
end
vj = vj';
quiver(pc(:,1),pc(:,2),vj(:,1),vj(:,2));


figure;
vr = [];
for i=1:length(vc)
    vr = [vr R*pc(i,:)'];
end
vr = vr';
quiver(pc(:,1),pc(:,2),vr(:,1),vr(:,2));


figure;
ve = [];
for i=1:length(vc)
    ve = [ve E*pc(i,:)'];
end
ve = ve';
quiver(pc(:,1),pc(:,2),ve(:,1),ve(:,2));


figure;
vd = [];
for i=1:length(vc)
    vd = [vd D*pc(i,:)'];
end
vd = vd';
quiver(pc(:,1),pc(:,2),vd(:,1),vd(:,2));


figure;
vs = [];
for i=1:length(vc)
    vs = [vs S*pc(i,:)'];
end
vs = vs';
quiver(pc(:,1),pc(:,2),vs(:,1),vs(:,2));


