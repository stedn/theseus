
bns = 51;
ll = 4;
rl = 40;

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


bp = '/Users/wmcfadden/activ_rec_sweep_b/';
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



if(size(allt)>0)
    save('actrec_meas','allt','allp','allg','alla','allf','alle','allc','alln')
end
