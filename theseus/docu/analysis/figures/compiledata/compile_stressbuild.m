
bns = 51;
ll = 10;
rl = 50;

allt = [];
allp = [];
allg = [];
allf = [];
alle = [];
allc = [];
allfe = [];
allfc = [];
alla = [];
allw = [];
alln = {};

bp = '/Users/wmcfadden/activ_llc/';
cd(bp);
files = dir;
files = {files.name};

for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            measureall
        end
    end
end

bp = '/Users/wmcfadden/activ_phi_sweep_a/';
cd(bp);
files = dir;
files = {files.name};

for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            measureall
        end
    end
end

bp = '/Users/wmcfadden/activ_phi_sweep_b/';
cd(bp);
files = dir;
files = {files.name};

for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            measureall
        end
    end
end


save('astress_meas','allt','allp','allg','alla','allf','alle','allc','allw','alln')
