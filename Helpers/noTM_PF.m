function noTM_PF(md)
%noTM_PF(md)
%
%   Runs CalculatePlacefields on the session, excluding treadmill runs. 
%
%   INPUT
%       md: Session entry
%
%   OUTPUT
%       See CalculatePlacefields. 
%

%% Get treadmill run indices. 
    cd(md.Location); 
    load('TimeCells.mat','TodayTreadmillLog'); 
    
    %Load data. 
    try
        load(fullfile(pwd,'Pos_align.mat'),'aviFrame');
    catch
        load('T2output.mat');
        [~,~,~,~,~,~,aviFrame] = AlignImagingToTracking(md.Pix2CM,FT,0);
    end
    
    %Treadmill runs. 
    inds = getTreadmillEpochs(TodayTreadmillLog,aviFrame); 
    nLaps = size(inds,1); 
    
    %Build index vector. 
    excludeme = []; 
    for i=1:nLaps
        excludeme = [excludeme, inds(i,1):inds(i,2)];
    end
    
%% Run CalculatePlacefields. 
    CalculatePlacefields(md,'exclude_frames',excludeme,'minspeed',3); 
end