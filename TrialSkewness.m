function skew = TrialSkewness(raster)
%
%
%   Rudimentary way to quantify when in the session a time cell comes
%   online. Takes the mean of trials active then divides by the total
%   number of trials. 

%%
    nLaps = size(raster,1);
    [trials,~] = find(raster); 
    
    meanTrial = mean(unique(trials)); 
    skew = meanTrial/nLaps;
end