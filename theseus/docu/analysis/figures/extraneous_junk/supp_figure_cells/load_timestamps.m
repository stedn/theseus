stamps = []
for i =1:200
    d=imfinfo(['/Volumes/will/stress/try3/fov1_' sprintf('%03d',i)  '.tif']);
    token = strsplit(d.DateTime, ':');
    stamps = [stamps;str2num(token{end})]; 
end
plot(stamps)