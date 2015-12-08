function TodayTreadmillLog = getTodayTreadmillLog(animal,date,session)
%TodayTreadmilllog = getTodayTreadmillLog(animal,date)
%
%   Gets the treadmill data and processes it for the given animal and time.
%   
%   WORKFLOW: MakeTreadmillMD -> getTodayTreadmillLog 
%
%   INPUTS
%       animal: String specifying the mouse (e.g., GCamp6f_45_treadmill). 
%
%       date: String specifying the date. (e.g., 11_20_2015).
%
%   OUTPUT
%       TodayTreadmillLog: Struct array with fields...
%           startts = vector for the timestamps associated with each
%           treadmill-on.
%
%           stopts = same for treadmill-off.
%
%           velocitysetting = setting for treadmill velocity.
%
%           delaysetting = setting for duration of delay. 
%
%           complete = logical for whether the full X seconds in
%           velocitysetting was run. 0 if the treadmill was stopped
%           manually.
%
%           speedEnforce = comes from TreadmillLog. 
%
%           numRuns = number of treadmill runs, complete and incomplete. 
%
%           StartTime = clock time of treadmill-on. 
%
%           StopTime = same for treadmill off. 
%
  
%% Find necessary indices. 
    LoadTreadmillMD;
    loadMD; 
    
    %Reference TreadmillMD. 
    animalInd = find(strcmp(animal,{TreadmillMD.Animal}));
    sessions = cell2mat({MD.Session});
    MDInd = find(strcmp(animal,{MD.Animal}) & strcmp(date,{MD.Date}) & sessions == session); 
    
    %Get file directory. 
    file = TreadmillMD(animalInd).File; 
    
    %Load treadmill history data file. 
    load(file); 
    
    %Get the dates in the treadmill history file. 
    dateformat(1:length(history.sessions)) = {'mm_dd_yyyy'};
    dates = cellfun(@datestr,{history.sessions.date},dateformat,'unif',0);
    
    %Find today's date. 
    dateInd = find(strcmp(date,dates));
    
%% Get the log. 
    %Timestamps are recorded from when you hit start on the TreadmillTask
    %GUI(?).
    firstTS = history.sessions(dateInd).treadmill(1,1);
    TodayTreadmillLog.startts = history.sessions(dateInd).treadmill(:,1)-firstTS;  %Normalize. 
    TodayTreadmillLog.stopts = history.sessions(dateInd).treadmill(:,2)-firstTS;
    TodayTreadmillLog.velocitysetting = history.sessions(dateInd).treadmill(:,3);
    TodayTreadmillLog.delaysetting = history.sessions(dateInd).treadmill(:,4);
    TodayTreadmillLog.complete = history.sessions(dateInd).treadmill(:,5);
    TodayTreadmillLog.RecordStartTime = MD(MDInd).RecordStartTime;
    
    %Get the clock times associated with treadmill on/off epochs. 
    timeformat = 'HH:MM:SS';    
    treadmilltimes = cell(length(history.sessions(dateInd).TreadmillLog(:,1)),1);
    for t=1:length(history.sessions(dateInd).TreadmillLog(:,1))
        treadmilltimes{t} = datestr(history.sessions(dateInd).TreadmillLog(t,1),timeformat);
    end
    speedEnforce = history.sessions(dateInd).TreadmillLog(:,2);
    
    TodayTreadmillLog.speedEnforce = speedEnforce(speedEnforce~=0);
    TodayTreadmillLog.numRuns = length(TodayTreadmillLog.complete); 
    
%% Get the clock time of the first treadmill-on. 
    %Duration of first run. 
    firstRunDuration = round(diff([TodayTreadmillLog.startts(1),TodayTreadmillLog.stopts(1)]));
    clockDurations = nan(length(treadmilltimes),3); 
    
    %Get all the treadmill durations. 
    for thisT=2:length(treadmilltimes)
        clockDurations(thisT-1,2) = etime(datevec(treadmilltimes{thisT}),datevec(treadmilltimes{thisT-1}));
        clockDurations(thisT-1,1) = clockDurations(thisT-1,2) - 1;
        clockDurations(thisT-1,3) = clockDurations(thisT-1,2) + 1;
    end
    
    %Index of first clock duration that matches that of the timestamp
    %durations. 
    [T,~] = find(ismember(clockDurations,firstRunDuration),1,'first');
    
    %First clock time of treadmill-on.
    TodayTreadmillLog.firstTreadmillOn = treadmilltimes{T};

end