function skew = TrialSkewness(raster,sigCurve)
%
%
%   Rudimentary way to quantify when in the session a time cell comes
%   online. Takes the mean of trials active then divides by the total
%   number of trials. 

%%
    nLaps = size(raster,1);
    [trials,~] = find(raster(:,sigCurve)>0); 
    
    meanTrial = mean(unique(trials)); 
    skew = meanTrial/nLaps;
end