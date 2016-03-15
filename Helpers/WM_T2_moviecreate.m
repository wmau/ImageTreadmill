% Script to take any number of movies you select and create all the imaging
% files necessary to run Tenaspis V2: DFMovie.h5 and SLPDF.h5

%% Filter specifications.
filter_type = 'circular';
filter_pixel_radius = 3;
LPfilter_pixel_radius = 20;

currDir = pwd; 

%% Step 0: Initialize.
mosaic.terminate;
mosaic.initialize;

%% Step 1: Select files to load and check for smoothed movies. 
%Motion corrected, cropped, un-smoothed movie inside the folder. 
MotCorrFiles = file_select_batch('*.mat'); 

nFiles = length(MotCorrFiles);      %Number of files. 
as3 = zeros(nFiles,1);              %Logical, 3-pixel smooth already done?
as20 = zeros(nFiles,1);             %Logical, 20-pixel smooth already done? 
amp = zeros(nFiles,1);              %Logical, minimum projection already done? 
adff = zeros(nFiles,1);             %Logical, DF/F already done? 

for i=1:nFiles
    %Session folder. 
    MotCorrFiles(i).sessionFolder = fileparts(MotCorrFiles(i).folder(1:end-1)); 
    
    %Is 3-pixel smooth already done? 
    as3(i) = exist(fullfile(MotCorrFiles(i).sessionFolder,'ICmovie_smoothed-Objects'),'dir') | ...
             exist(fullfile(MotCorrFiles(i).sessionFolder,'ICmovie_smooth_circular_3-Objects'),'dir'); 
         
    %Is 20-pixel smooth already done? 
    as20(i) = exist(fullfile(MotCorrFiles(i).sessionFolder,'LPmovie_circular_20-Objects'),'dir'); 
    
    %Is minimum projection already done? 
    amp(i) = exist(fullfile(MotCorrFiles(i).sessionFolder,'ICmovie_min_proj.tif'),'file'); 
    
    %Is DF/F already done? 
    adff(i) = exist(fullfile(MotCorrFiles(i).sessionFolder,'DFF-Objects'),'dir'); 
end

%% Start processing movies. 
for i=1:nFiles
    MotCorrMovie = mosaic.loadObjects(MotCorrFiles(i).path);
    cd(MotCorrFiles(i).sessionFolder); 
    
    %Directories for saving. 
    threePixSmooth.name = fullfile(pwd,['ICmovie_smooth_' filter_type '_' num2str(filter_pixel_radius)]);  
    threePixSmooth.folder = [threePixSmooth.name,'-Objects'];
    threePixSmooth.mat = [threePixSmooth.name,'.mat'];
    
    DFF.name = fullfile(pwd,'DFF'); 
    DFF.folder = [DFF.name,'-Objects'];
    DFF.mat = [DFF.name,'.mat'];
      
    LP.name = fullfile(pwd,['LPmovie_' filter_type '_' num2str(LPfilter_pixel_radius)]); 
    LP.folder = [LP.name,'-Objects'];
    LP.mat = [LP.name,'.mat'];
    
%% Minimum projection. 
    if ~amp(i)  %Do min project
        ICmovie_min_proj = mosaic.projectMovie(MotCorrMovie,'projectionType','Minimum'); 
        mosaic.saveImageTiff(ICmovie_min_proj,'ICmovie_min_proj.tif');
        
    else        %Load min project
        disp('Minimum projection already done! Delete ICmovie_min_proj.tif to rerun.'); 
        ICmovie_min_proj = mosaic.loadImage(fullfile(pwd,'ICmovie_min_proj.tif'));
    end
    
%% 3-pixel radius disc filter. 
    if ~as3(i) 
        disp(['Performing ' num2str(filter_pixel_radius) ' pixel disc smoothing of motion corrected movie']);
        %Perform filter.
        threePixSmooth.movie = mosaic.filterMovie(inputMovie,'filterType', filter_type,...
            'filterSize',filter_pixel_radius*2); 
        %Save. 
        mosaic.saveOneObject(threePixSmooth.movie,threePixSmooth.name);    
        
    else
        disp([num2str(filter_pixel_radius),' pixel smooth already done.']); 
        
        %Load. 
        if ~adff(i)
            threePixSmooth.movie = mosaic.loadObjects(threePixSmooth.mat); 
        end
    end
    
    %Folder contents. 
    cd(threePixSmooth.folder); 
    threePixSmooth.h5 = fullfile(pwd,ls('*.h5'));
    threePixSmooth.foldermat = fullfile(pwd,ls('*.mat')); 
    cd(MotCorrFiles(i).sessionFolder); 
    
%% DF/F
    if ~adff(i)
        disp('Creating DFF movie');
        %Do DF/F.
        DFF.movie = mosaic.normalizeMovie(threePixSmooth.movie,...
            'method', '(f-f0)/f0','f0Image',ICmovie_min_proj);
        %Save.
        mosaic.saveOneObject(DFF.movie,'DFFmovie'); 
    else
        disp('DFF movie already ran.'); 
    end
    
    clear DFF threePixSmoothMovie;
    
%% 20-pixel radius disc filter. 
    if ~as20(i)
        disp(['Performing ' num2str(LPfilter_pixel_radius) ' pixel disc smoothing of motion corrected movie for LPmovie']);
        %Perform filter.
        LP.movie = mosaic.filterMovie(MotCorrMovie,'filterType', filter_type,...
            'filterSize',LPfilter_pixel_radius*2);         
        %Save.
        mosaic.saveOneObject(LP.movie,LP.name); 
    else
        disp([num2str(LPfilter_pixel_radius), ' pixel smooth already done.']); 
    end
    
    %Folder contents. 
    cd(LP.folder); 
    LP.h5 = fullfile(pwd,ls('*.h5'));
    LP.foldermat = fullfile(pwd,ls('*.mat')); 
    cd(MotCorrFiles(i).sessionFolder); 
    
%% TS_Lowpass_Divide
    disp('Creating TS Lowpass Divide movie')
    TS_Lowpass_Divide(threePixSmooth.h5,LP.h5);
    
    % Cleanup everything
    disp(['Deleting smoothed movie 1: ',threePixSmooth.h5]);
    delete(threePixSmooth.h5);
    delete(threePixSmooth.foldermat);
    delete(threePixSmooth.mat); 
    
    disp(['Deleting smoothed movie 2: ',LP.h5]);
    delete(LP.h5);
    delete(LP.foldermat);
    delete(LP.mat); 
    
end