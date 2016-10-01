function [r,lags]=xcorr_by_laps(triggerRaster,targetRaster)
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
%       lags: Vector of lags. 
%

%% Perform cross correlation.
    nLaps = size(triggerRaster,1); 
    
    r = zeros(nLaps,399); 
    for l=1:nLaps
        [r(l,:),lags] = xcorr(triggerRaster(l,:),targetRaster(l,:),'coeff'); 
    end
    
    %Lags is in frames. 
end