function [r,lags]=xcorr_by_laps(triggerRaster,targetRaster)
%
%
%

%%
    nLaps = size(triggerRaster,1); 
    
    r = zeros(nLaps,399); 
    for l=1:nLaps
        [r(l,:),lags] = xcorr(triggerRaster(l,:),targetRaster(l,:),'coeff'); 
    end
    
    %Lags is in frames. 
end