function tuningStatus = TCRemap(ref,ssn)
%tuningStatus = TCRemap(ref,ssn)
%
%   Finds how time cells in one sessions are consistent across other days. 
%
%   INPUTS    
%       ref: MD entry, the session containing the time cells that you want
%       to track remapping across days. 
%
%       ssn: MD entry, other sessions that you are comparing to.

%   OUTPUT
%       sametuning: NxS matrix (N=neurons, S=sessions). Values can be: 
%           NaN: First column and neurons that were not mapped. 
%           -1: No longer time cells. 
%           0: Remapped to a different time. 
%           1: Same response curve. 
%      

%% Compile time cell data and get mapping. 
    %Partition the session data. 
    ssns = [ref, ssn];              %First session of MD is the base session.
    nSessions = length(ssns);       %Number of sessions, including base and comparison sessions.
 
    %Get all the time cell data from each session. 
    DATA = CompileMultiSessionData(ssns,{'timecells','curves','t'});
    TIMECELLS = DATA.timecells; 
    CURVES = DATA.curves; 
    Ts = cell2mat(DATA.t);
    nNeurons = length(DATA.curves{1}.tuning);    
    
%% Find the indices in MAP corresponding to sessions of interest.
    mapMD = getMapMD(ref);
    matchMat = msMatchCells(mapMD,ssns,TIMECELLS{1},true);  
    
%% Indicate whether a neuron has remapped or not. 
    %Find time resolution of the tuning curves. 
    tResolution = zeros(nSessions,1);
    for s=1:nSessions
        nBins = length(CURVES{s}.tuning{1});
        tResolution(s) = Ts(s)/nBins;
    end
    
    %If there are multiple time resolutions, throw an error.
    assert(length(unique(tResolution))==1,...
        'Different time resolutions found in tuning curves! Rerun FindTimeCells!');
    tResolution = unique(tResolution);              %seconds. 
    
    window = .75;                                   %seconds.
    binWindow = round(1/tResolution*window);        %bins.
    tuningStatus = nan(nNeurons,nSessions);         %Preallocate. 
    nTCs = size(matchMat,1);
    
    %For each time cell in the base session...
    for i = 1:nTCs
        n1 = matchMat(i,1);
        baseSig = CURVES{1}.sig{n1};        %Significance curve for base session.

        %Find humps in the significance curve (essentially translates to
        %the tuning curve).
        ccBase = bwconncomp(baseSig);               %Using connected components. 
        nHumps = length(ccBase.PixelIdxList);       %Number of humps in the curve. 

        %Get bin numbers corresponding to the middle of the hump. 
        baseBins = zeros(nHumps,1); 
        for h=1:nHumps
            baseBins(h) = round(mean(ccBase.PixelIdxList{h}));
        end

        %For each comparison session...
        for s=2:nSessions
            %Second neuron.
            n2 = matchMat(i,s);
            
            %Significance curve of comparison session.
            compSig = CURVES{s}.sig{n2};

            %Find humps in the significance curve of the comparison
            %session. 
            ccComp = bwconncomp(compSig);
            nHumps = length(ccComp.PixelIdxList); 

            %Get bin numbers corresponding to the middle of the hump. 
            compBins = zeros(nHumps,1);
            for h=1:nHumps
                compBins(h) = round(mean(ccComp.PixelIdxList{h}));
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
            if isempty(b) || ~ismember(n2,TIMECELLS{s})
                tuningStatus(n1,s) = -1; 
            %If sigcurve of the two sessions are within binWindow
            %bins, 1.               
            elseif any(sameTimeField)
                tuningStatus(n1,s) = 1; 
            %If comparison sigcurve is still a time cell, but encodes a
            %different time, 0.
            elseif any(differentTimeField)
                tuningStatus(n1,s) = 0;                                 
            end

        end
        
    end
    
    
end