function inds = getTreadmillEpochs(TodayTreadmillLog,time_interp)
% inds = getTreadmillEpochs(TodayTreadmillLog,time_interp)
%
%   Gets the indices that will reference FT for when the start and stop
%   times of treadmill runs occur. 
%
%   INPUTS
%       TodayTreadmillLog: output from getTodayTreadmillLog.
%
%       time_interp: aligned timestamp data from Pos.mat. 
%
%   OUTPUT
%       inds: Lx2 matrix (L=number of laps) where its values are the
%       indices for FT of the start (first column) and end (second column)
%       of treadmill runs. 
%

%% Find closest timestamp for treadmill runs. 
    inds = zeros(TodayTreadmillLog.numRuns,2);
    for thisRun = 1:TodayTreadmillLog.numRuns
        inds(thisRun,1) = findclosest(TodayTreadmillLog.startts(thisRun),time_interp);
        inds(thisRun,2) = findclosest(TodayTreadmillLog.stopts(thisRun),time_interp);
    end
    
end