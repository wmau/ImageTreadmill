function shuffled = permuteTime(raster)
%shuffled = permuteTime(raster)
%
%
    
%%
    [nLaps,nBins] = size(raster);
    r = rem(randi([1,200],nLaps,1),nBins);
    wrapper = [raster,raster]; 
    shuffled = wrapper(bsxfun(@plus,bsxfun(@plus,nBins-r,0:nBins-1)*nLaps,(1:nLaps)'));
end