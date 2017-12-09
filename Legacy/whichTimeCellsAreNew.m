function [remapped,MAP,MAPcols,TIMECELLS,CURVES] = whichTimeCellsAreNew(tday,yday,MAPlocation,Ts)
%[newTCs,MAP,MAPcols,TIMECELLS,CURVES] = whichTimeCellsAreNew(tday,yday,MAPlocation,Ts)
%   
%   Take two recording sessions reasonably close in time (say one day and
%   the day before). Across these two sessions there will be neurons that
%   remap (while still remaining as time cells) and new neurons that start
%   to encode time. This function looks for all such neurons and outputs an
%   index for batch_session_map.map for those neurons. 
%
%   INPUTS
%       tday: MD entry for 'base' recording session. Effectively this means
%       that the output will be neurons from this session that gained a new
%       time field. 
%
%       yday: MD entry for 'comparison' recording session from which the
%       change in temporal response will be measured. 
%
%       MAPlocation: Directory with batch_session_map.
%
%       Ts: Vector of delay durations, two elements in this function. 
%
%   OUTPUTS
%       newTCs: Vector of indices for MAP, neurons that gained a time
%       response in tday. 
%
%       MAP: Matrix containing neuron mappings. 
%
%       MAPcols: Column indices for MAP corresponding to the sessions
%       indicated in the inputs. 
%
%       TIMECELLS: Two cells containing time cells referencing the
%       respective FTs of that session. First cell array is the yday
%       TimeCells. 
%
%       CURVES: Two cell containing curves from FindTimeCells. 
%
 
%% Get MAP and session columns. 
    %Let's organize in chronological order. 
    dates = {yday.Date, tday.Date}; 
    sessionNums = [yday.Session, tday.Session];

    load(fullfile(MAPlocation.Location,'batch_session_map.mat')); 
    
    [MAP,MAPcols] = FilterMAPDates(batch_session_map,dates,sessionNums);
    
%% Get new time response neurons. 
    %Find neurons from yesterday that have different time tuning today (0
    %in the second column of sametuning). 
    [sametuning,TIMECELLS,CURVES] = TimeCellRemapRate(MAPlocation,yday,tday,Ts);
    
    %Neurons that remap. 
    remap = find(sametuning(:,2)==0);                               %Indices, referencing TimeCells{1}. 
    newTuning = find(ismember(MAP(:,MAPcols(1)),TIMECELLS{1}(remap))); %Indices, referencing MAP. 
    
    %Find time cells from tday and yday. 
    ydayTC = find(ismember(MAP(:,MAPcols(1)),TIMECELLS{1}));        %Indices, referencing MAP.
    tdayTC = find(ismember(MAP(:,MAPcols(2)),TIMECELLS{2}));        %Indices, referencing MAP.
    
    %Both new time cells and time cells that have remapped. 
    newTCs = tdayTC(~ismember(tdayTC,ydayTC)); 
    
    remapped.new = newTCs; 
    remapped.tuning = newTuning;
    
end