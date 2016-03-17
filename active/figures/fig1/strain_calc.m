bp = '/Users/wmcfadden/Desktop/jons data/examples/sym/';
code = 'SYM_3_TRACK_retry.txt';
in = importdata([bp code],'\t',1);
data = in.data;

tstep = 2;

p = data(data(:,3)==1,4:5);
stoD = [0];
stot = [0];
for t = 2:max(data(:,3))
    subind = data(:,3)==t;
    pt = data(subind,4:5);
    v = pt-p;
    [p0,v0,J,R,E,D,S] = strainrate(p,v);
    stoD = [stoD; D];
    stot = [stot; (t-1)*tstep];
end

% plot(stot,stoD)
cftool(stot,1+stoD)