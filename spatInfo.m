function spatInfo(md)
%
%
%

%% Load.
    cd(md.Location);
    load('PlaceMaps.mat','TMap_unsmoothed','OccMap','isrunning');
    load('Pos_align.mat','FT');
    
    nNeurons = length(TMap_unsmoothed); 
    nFrames = sum(isrunning);
    flatOccMap = OccMap(:);
    p = flatOccMap ./ nFrames;
    thresh = 4/nFrames;
    bad = p < thresh;
    p(bad) = [];
    lambda = sum(FT(:,isrunning))./nFrames;
    
    spatialI = zeros(1,nNeurons);
    for n=1:nNeurons
        lambda_i = TMap_unsmoothed{n}(:); 
        lambda_i(bad) = [];
        
        spatialI(n) = nansum(p.*lambda_i.*log2(lambda_i./lambda(n)));
    end
    
    save('SpatialInfo.mat','spatialI');
end