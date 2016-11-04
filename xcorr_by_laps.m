function [r,lags]=xcorr_by_laps(triggerRaster,targetRaster,corrType)
%[r,lags]=xcorr_by_laps(triggerRaster,targetRaster)
%
%   Performs lap by lap cross correlation of two rasters. 
%
%   INPUTS
%       triggerRaster and targetRaster: Logical arrays indicating onset
%       times for spiking during treadmlil runs. 
%
%   OUTPUTS
%       r: LxT (L=# of laps, T=# of lag bins) matrix for cross
%       correlations at each lap.
%
%       lags: Vector of lags, in frames. 
%

%% Perform cross correlation.
    nLaps = size(triggerRaster,1); 
    maxlag_seconds = 4;             %Default hard code.
    maxlag = maxlag_seconds*20;     %Translate into frames. 
    
    %XCorr. 
    r = zeros(nLaps,maxlag*2+1); 
    for l=1:nLaps
        [r(l,:),lags] = xcorr(triggerRaster(l,:),targetRaster(l,:),maxlag,corrType); 
    end
    
    %Lags is in frames. 
end