function combineIndexCSVwithFakeDVT(videoFile)
%
%
%

%%
    [~,base_name] = fileparts(videoFile);
    vidobj = VideoReader(videoFile);
    num_frames = vidobj.Duration*vidobj.FrameRate;
    fakeDVT = zeros(num_frames, 5);
    fakeDVT(1:num_frames,1) = 1:num_frames;
    
    csvFile = csvread([base_name,'.csv']); 
    fakeDVT(1:num_frames,2) = csvFile(2:end,2);
    
    csvwrite([base_name '.DVT'],fakeDVT);
end