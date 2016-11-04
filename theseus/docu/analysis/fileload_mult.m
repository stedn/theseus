bp = '/Users/wmcfadden/activ_stress_sweep/';
cd(bp);
files = dir;
files = {files.name};
for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            fileload_one
        end
    end
end
