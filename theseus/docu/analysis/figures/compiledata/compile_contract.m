bp = '/Users/wmcfadden/activ_free_sweep_nonl/';
cd(bp);
files = dir;
files = {files.name};

bns = 51;
ll = 10;
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
for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            measureall
        end
    end
end
save('contract_meas','allt','allp','allg','alla','allf','alle','allc','allw','alln')
