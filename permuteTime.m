function shuffled = permuteTime(raster)
%shuffled = permuteTime(raster)
%
%   Performs a randomized card shuffle permutation for each row in raster. 
%
%   INPUT
%       raster: Trials x Time matrix of anything.
%
%   OUTPUT
%       shuffled: Same matrix, after shuffling elements in each row by a
%       random value. 
%
    
%% Shuffle time. 
    [nLaps,nBins] = size(raster);
    r = rem(randi([1,200],nLaps,1),nBins);
    wrapper = [raster,raster]; 
    shuffled = wrapper(bsxfun(@plus,bsxfun(@plus,nBins-r,0:nBins-1)*nLaps,(1:nLaps)'));
end