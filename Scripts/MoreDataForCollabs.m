directory = 'B:\Treadmill project'; 
files = dir(directory); 

cd(directory);
c = 1; 
for i=3:6
    cd(files(i).name);
    
    date_folders = dir;
    
    for j=3:6
        cd(date_folders(j).name);
        
        folder = dir;
        cd(folder(3).name); 
        
        load('Pos_align.mat', 'RawTrace');
        load('TreadmillTraces.mat', 'RawTrdmll'); 
        
        [x,Date] = fileparts(folder(3).folder);
        [x,Animal] = fileparts(x);
        
        DATA(c).Date = Date;
        DATA(c).Animal = Animal;
        DATA(c).TreadmillTraces = RawTrdmll;
        DATA(c).Traces = RawTrace; 
        
        c = c+1;
        cd ..
        cd ..
    end
    
    cd ..
end