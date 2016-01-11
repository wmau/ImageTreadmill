a=dir;

for i=11:20
    try
        cd(fullfile(a(i).name,'1 - Homecage'));
    catch
        cd(fullfile(a(i).name,'2 - Homecage'));
    end
    XMLfile = dir('*.xml');
    WM_batch_preprocess(fullfile(pwd,XMLfile.name));
    cd ..
    cd ..
end