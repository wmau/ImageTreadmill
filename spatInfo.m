function spatInfo(md)
%
%
%

%% Load.
    cd(md.Location);
    load('PlaceMaps.mat','TMap_unsmoothed','RunOccMap','isrunning','cmperbin',...
        'frames_use_ind');
    load('Pos_align.mat','FT');
    good = frames_use_ind & isrunning;
    nFrames = sum(good); 
    
%% Smooth occupancy map.
    p = RunOccMap./nFrames;
    p = p(:);
    
    nNeurons = length(TMap_unsmoothed); 
    thresh = 5; %frames
    bad = RunOccMap < 4;
    p(bad) = [];
    lambda = sum(FT(:,good),2)'./nFrames;
    
    spatialI = zeros(1,nNeurons);
    mutInfo = zeros(1,nNeurons);
    for n=1:nNeurons
        lambda_i = TMap_unsmoothed{n}(:); 
        lambda_i(bad) = [];

        spatialI(n) = nansum(p.*lambda_i.*log2(lambda_i./lambda(n)));
        
        mutInfo(n) = calc_mutual_information(TMap_unsmoothed{n},RunOccMap);
    end
    
    spec = spatialI./lambda;
    
    save('SpatialInfo.mat','spec','spatialI','lambda','mutInfo');
end