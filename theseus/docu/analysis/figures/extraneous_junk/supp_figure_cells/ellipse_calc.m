bp = '/Users/wmcfadden/Desktop/stressed_out/sequences/';
codes = {'wt1_1.txt' 'wt1_2.txt' 'wt1_3.txt' 'wt3_1.txt' 'wt4_1.txt' 'wt5_1.txt' 'wt5_2.txt'};
tsteps = [20 20 10 20 8 3 3];
tau = [];
for i = 1:length(codes)
    code = codes{i}
    in = importdata([bp code],'\t',1);
    data = in.data;

    tstep = tsteps(i);


    stoa = [];
    stob = [];
    stot = [];
    for t = 1:max(data(:,3))
        subind = data(:,3)==t;
        pt = data(subind,4:5);
        ot = fit_ellipse(pt(:,1),pt(:,2));
        a = max(ot.a,ot.b);
        b = min(ot.a,ot.b);
        stoa = [stoa; a];
        stob = [stob; b];
        stot = [stot; (t-1)*tstep];
    end
    r = stoa./stob;
    f=fit(stot,r,'a*exp(-b*x)+1','StartPoint',[1 0.01]);
    tau = [tau; 1/f.b];
    plot(stot,r);
    hold on
    plot(f);
end
% plot(stot,stoD)
