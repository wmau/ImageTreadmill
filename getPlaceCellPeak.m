function peaks = getPlaceCellPeak(md)
%peaks = getPlaceCellPeak(md)
%
%   Get peaks for place cells in linearized space from 2.5 cm bins. 
%
%   INPUT
%       md: session entry.
%   
%   OUTPUT
%       peaks: N-length vector containing peaks of tuning curve. 
%
    
%% Get peaks.
    cd(md.Location);
    load('SpatialTraces.mat','curve'); 
        
    [~,peaks] = max(curve,[],2); 
    peaks = peaks - 22;                         %Bins 1-22 contain the treadmill.

    peaks(peaks<=0) = peaks(peaks<=0) + 80;     %Basically a circular permutation around the maze. 
end