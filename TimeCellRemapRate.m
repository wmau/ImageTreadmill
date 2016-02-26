function [sametuning,TIMECELLS,CURVES,MAP,MAPcols] = TimeCellRemapRate(mapMD,base,comp,Ts)
%[sametuning,TIMECELLS,CURVES] =
%TimeCellRemapRate(batch_session_map,base,comp,Ts)
%
%   Finds how time cells in one sessions are consistent across other days. 
%
%   INPUTS
%       batch_session_map: Output from neuron_reg_batch. 
%
%       base: MD entry, the session containing the time cells that you want
%       to track remapping across days. 
%
%       comp: MD entry, other sessions that you are comparing to.
%
%       Ts: Vector, same length as base+comp, containing delay durations
%       for each session.
%
%   OUTPUTS
%       sametuning: NxS matrix (N=neurons, S=sessions). Values can be: 
%           NaN: First column and neurons that were not mapped. 
%           -1: No longer time cells. 
%           0: Remapped to a different time. 
%           1: Same response curve. 
%       
%       TIMECELLS: Sx1 cell containing vectors of time cell indices. 
%
%       CURVES: Sx1 cell containing curves from FindTimeCells. 
%

%% Compile time cell data and get mapping. 
    %Partition the session data. 
    MD = [base, comp];              %First session of MD is the base session.
    nSessions = length(MD);         %Number of sessions, including base and comparison sessions.
    dates = {MD.Date};              %All dates. 
    sessions = [MD.Session];        %Session numbers. 
      
    %Get all the time cell data from each session. 
    [TIMECELLS,~,CURVES,~,~] = CompileTimeCellData(MD,Ts);
    nTimeCells = length(TIMECELLS{1}); 
    
    %Get neuron mapping. 
    load(fullfile(mapMD.Location,'batch_session_map.mat')); 
    
%% Find the indices in MAP corresponding to sessions of interest.
    regDates = {batch_session_map.session.Date};
    regSessions = [batch_session_map.session.Session]; 
    
    [MAP,MAPcols] = FilterMAPDates(batch_session_map,dates,sessions);
   
    %Rows on MAP corresponding to where time cells are found in the base
    %session. 
    MAProws = zeros(nTimeCells,1);
    for i=1:nTimeCells
        row = find(ismember(MAP(:,MAPcols(1)),TIMECELLS{1}(i)));
        if ~isempty(row), MAProws(i) = row; end
    end
    
%% Indicate whether a neuron has remapped or not. 
    %Find time resolution of the tuning curves. 
    tResolution = zeros(nSessions,1);
    for thisSession=1:nSessions
        nBins = length(CURVES{thisSession}.tuning{1});
        tResolution(thisSession) = Ts(thisSession)/nBins;
    end
    
    %If there are multiple time resolutions, throw an error.
    assert(length(unique(tResolution))==1,...
        'Different time resolutions found in tuning curves! Rerun FindTimeCells!');
    tResolution = unique(tResolution);              %seconds. 
    
    window = 0.5;                                   %seconds.
    binWindow = round(1/tResolution*window);        %bins.
    sametuning = nan(length(MAProws),nSessions);    %Preallocate. 
    ind = 1;
    
    %For each time cell in the base session...
    for thisRow=MAProws'
        if thisRow
            neurons = MAP(thisRow,MAPcols);             %Vector indexing MAP. 
            neurons(isnan(neurons)) = 0;                %Get rid of NaNs. 
            baseSig = CURVES{1}.sig{neurons(1)};        %Significance curve for base session.

            %Find humps in the significance curve (essentially translates to
            %the tuning curve).
            ccBase = bwconncomp(baseSig);               %Using connected components. 
            nHumps = length(ccBase.PixelIdxList);       %Number of humps in the curve. 

            %Get bin numbers corresponding to the middle of the hump. .
            baseBins = zeros(nHumps,1); 
            for i=1:nHumps
                baseBins(i) = round(mean(ccBase.PixelIdxList{i}));
            end

            %For each comparison session...
            for thisSession=2:nSessions
                if neurons(thisSession)                 %If the neuron was detected in neuron_reg_batch.
                    %Significance curve of comparison session.
                    compSig = CURVES{thisSession}.sig{neurons(thisSession)};

                    %Find humps in the significance curve of the comparison
                    %session. 
                    ccComp = bwconncomp(compSig);
                    nHumps = length(ccComp.PixelIdxList); 

                    %Get bin numbers corresponding to the middle of the hump. 
                    compBins = zeros(nHumps,1);
                    for i=1:nHumps
                        compBins(i) = round(mean(ccComp.PixelIdxList{i}));
                    end

                    %All possible combinations. 
                    [b,c] = meshgrid(baseBins,compBins); 

                    %Grid of time field comparisons within the binWindow. 
                    sameTimeField = abs(b-c) <= binWindow; 
                    sameTimeField = sameTimeField(:);               %This needs to be vectorized for any() to work. 
                    differentTimeField = abs(b-c) > binWindow;
                    differentTimeField = differentTimeField(:);              

                    %Subtract. Order of these if statements matters here! 
                    %If no longer a time cell from lack of a response or not
                    %listed in TIMECELLS, -1. 
                    if isempty(b) || ~ismember(neurons(thisSession),TIMECELLS{thisSession})
                        sametuning(ind,thisSession) = -1; 
                    %If sigcurve of the two sessions are within binWindow
                    %bins, 1.               
                    elseif any(sameTimeField)
                        sametuning(ind,thisSession) = 1; 
                    %If comparison sigcurve is still a time cell, but encodes a
                    %different time, 0.
                    elseif any(differentTimeField)
                        sametuning(ind,thisSession) = 0;                                 
                    end

                end

            end
        end
        
        ind = ind+1;
    end
    
    
end