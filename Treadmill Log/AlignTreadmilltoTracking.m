function TodayTreadmillLog = AlignTreadmilltoTracking(TodayTreadmillLog,RecordStartTime)
%TodayTreadmillLog = AlignTreadmilltoTracking(TodayTreadmillLog,RecordStartTime)
%
%   WORKFLOW: MakeTreadmillMD -> getTodayTreadmillLog ->
%   AlignTreadmilltoTracking
%   
%   Align the treadmill on/off timestamps to the tracking timestamps by
%   finding the difference in seconds between the first treadmill-on and
%   the time recording starts. Then add those to the treadmill timestamps
%   after having made the first timestamp 0 upstream. 
%
%   INPUTS
%       TodayTreadmillLog: output from getTodayTreadmillLog
%
%       RecordStartTime: clock time at which you started the tracking
%       recording. 
%
%   OUTPUT
%       TodayTreadmillLog: same as input, but with timestamps aligned to
%       tracking.
%

%% Alignment. 
    %Determine the offset in seconds of the recording start time and the
    %first treadmill run. Then add this number to all timestamps in
    %TodayTreadmillLog. 
    offset = etime(datevec(TodayTreadmillLog.firstTreadmillOn),datevec(RecordStartTime));
    
    TodayTreadmillLog.startts = TodayTreadmillLog.startts + offset; 
    TodayTreadmillLog.stopts = TodayTreadmillLog.stopts + offset; 
    
end