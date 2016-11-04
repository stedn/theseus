
bns = 51;
ll = 10;
rl = 50;

load('../../../data/astress_meas.mat')

bp = '/Users/wmcfadden/activ_stress_sweep/';
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


save('astress_meas_correct','allt','allp','allg','alla','allf','alle','allc','allw','alln')
