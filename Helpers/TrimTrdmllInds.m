function [inds,nRuns] = TrimTrdmllInds(TodayTreadmillLog,T)
%[inds,nRuns] = TrimTrdmllInds(TodayTreadmillLog,T)
%
%   Makes the treadmill indices a fixed length (20*T - 1). 
%
%   INPUTS
%       TodayTreadmillLog in TimeCells.mat.
%
%       T: Run duration in seconds.
%
%   OUTPUTS
%       inds: TodayTreadmillLog.inds, equal lengths.
%
%       nRuns: Number of complete runs. 
%

%% Trim indices. 
    inds = TodayTreadmillLog.inds;
    inds = inds(find(TodayTreadmillLog.complete),:);
    inds(:,2) = inds(:,1) + 20*T-1;
    
    nRuns = sum(TodayTreadmillLog.complete); 
    
end