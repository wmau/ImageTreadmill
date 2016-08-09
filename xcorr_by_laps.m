function xcorr_by_laps(triggerRaster,targetRaster,critLaps)
%
%
%

%%
    %Turn into numerics. 
    triggerRaster = single(triggerRaster); 
    targetRaster = single(targetRaster); 
    
    [triggerLaps,~] = find(triggerRaster); 
    [targetLaps,~] = find(targetRaster); 
    
    LapsBothActive = intersect(triggerLaps,targetLaps)';
    nLapsBothActive = length(LapsBothActive); 
    
    lapnum = 1; 
    r = zeros(nLapsBothActive,399); 
    if length(LapsBothActive) > critLaps
        for l=LapsBothActive
            [r(lapnum,:),lags] = xcorr(triggerRaster(l,:),targetRaster(l,:)); 
            
            plot(lags,r(lapnum,:)); hold on;
            
            lapnum = lapnum+1; 
        end
    end
    
   
end
            