function decodeError = TimeDecodeError(decodedTime)
%
%
%

%%
    [nBins,nRuns] = size(decodedTime); 
    t = linspace(0,10,nBins)';
    
    realTime = repmat(t,1,nRuns);
    decodeError = abs(decodedTime - realTime); 
end