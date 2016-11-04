function [b,n,s]=bindata_line(x,y,gx)
    % [b,n,s]=bindata(x,y,gx)
    % Bins y(x) onto b(gx), gx defining centers of the bins. NaNs ignored.
    % Optional return parameters are:
    % n: number of points in each bin   
    % s: standard deviation of data in each bin 
    % A.S.
    x1 = x(:,1);
    x2 = x(:,3);
    n=[];
    b=[];
    s=[];
    for g=gx
        subs = (x1<g&x2>g)|(x2<g&x1>g);
        n=[n sum(subs)];
        b=[b mean(y(subs))];
        s=[s std(y(subs))];
    end

end