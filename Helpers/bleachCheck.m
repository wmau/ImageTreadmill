function [m,sem,chunkID] = bleachCheck(ts,nBins)
%[m,sem,chunkID] = bleachCheck(ts,nBins)
%
%   Bins a fluorescence time series 

%%
    T = length(ts);
    chunkSize = round(T/nBins);
    chunkLims = 0:chunkSize:T;
    chunkLims(end) = T;
    
    [m,sem] = deal(zeros(1,nBins));
    chunkID = zeros(size(ts));
    for b=1:nBins
        inds = chunkLims(b)+1:chunkLims(b+1);
        
        m(b) = mean(ts(inds));
        sem(b) = std(ts(inds))./sqrt(length(inds));
        
        chunkID(inds) = b;
    end
end