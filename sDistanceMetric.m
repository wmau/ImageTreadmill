function [stability,mds] = sDistanceMetric(sortedPeaks)
%[stability,mds] = sDistanceMetric(sortedPeaks,mds,position)
%   
%   Computes the distance that the peak time tuning curve drifted since the
%   day indicated by the input argument position. 
%
%   INPUTS
%       sortedPeaks: This is from multiPastalkovaPlot. The first column
%       consists of the tuning curve peaks (in seconds). Subsequent columns
%       are the peaks for the same neurons on the sessions in the same
%       order as presented in the input argument for multiPastalkova plot,
%       compMDs.
%
%       position: Specifies which session 

%% Setup.
    %Get important variables and preallocate. 
    [nTCs,nSessions] = size(sortedPeaks); 
    stability = nan(nTCs,nSessions); 
  
    %Get the distance metric by subtracting the peaks of the base session
    %from every other session.
    stability(:,1) = zeros(nTCs,1); 
    for s=2:nSessions
        stability(:,s) = abs(sortedPeaks(:,s) - sortedPeaks(:,1)); 
    end
    
end