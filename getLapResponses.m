function [ratebylap,delays,x,y,time_interp] = getLapResponses(animal,date,sessionNum,FT,TodayTreadmillLog)
%ratebylap = getLapResponses(animal,date,sessionNum,FT,TodayTreadmillLog)
%
%   Get the lap by lap responses for each neuron during treadmill run. 
%
%   INPUTS
%       animal: animal ID (e.g., GCamp6f_45_treadmill).
%   
%       date: session date string (e.g., 11_20_2015).
%
%       sessionNum: session number. 
%
%       FT: from ProcOut. Not aligned.
%
%       TodayTreadmillLog: from getTodayTreadmillLog.
%
%   OUTPUT
%       ratebylap: LxBxN matrix (L=number of runs, B=number of bins, N =
%       number of neurons) containing histograms for every neuron time
%       responses on the treadmill. 
%
%       delays: Vector of length L containing the delay setting for each
%       lap. 
%
%       X&Y: Aligned position data. 
%       
%       time_interp: Interpolated and aligned timestamps. 
%

%% Preliminary stuff. 
    ChangeDirectory(animal,date,sessionNum);
    
    %Align FT. 
    [x,y,~,FT,~,~,~,time_interp] = AlignImagingToTracking(0.15,FT);
    
    %Get treadmill run epochs. 
    inds = getTreadmillEpochs(TodayTreadmillLog,time_interp);
    nFramesBetween = diff(inds,1,2);                        %Number of frames epochs span. 
    
%% Bin time responses. 
    %Initialize.  
    [nNeurons,~] = size(FT); 
    tResolution = 0.20;                                     %seconds
    nBins = TodayTreadmillLog.delaysetting/tResolution;     %vector specifying number of bins per lap
    nComplete = sum(TodayTreadmillLog.complete); 
    completeLaps = find(TodayTreadmillLog.complete); 
    ratebylap = nan(nComplete,max(nBins),nNeurons); 
    
    p = ProgressBar(nNeurons);
    for thisNeuron=1:nNeurons
        for thisLap=1:nComplete
            lapnum = completeLaps(thisLap);
            tStart = 0;
            tEnd = TodayTreadmillLog.stopts(lapnum) - TodayTreadmillLog.startts(lapnum); 
            %Time vector. 
            t = linspace(tStart,tEnd,nFramesBetween(lapnum)+1); 

            %Times where there was a spike. 
            tspk = t(FT(thisNeuron,inds(lapnum,1):inds(lapnum,2))==1); 

            %Edges for histogram.
            edges = linspace(0,TodayTreadmillLog.delaysetting(lapnum),nBins(lapnum));

            %Make histogram.
            [ratebylap(thisLap,1:nBins(lapnum),thisNeuron),~] = hist(tspk,edges); 
        end
        p.progress;
    end
    p.stop;

    delays = TodayTreadmillLog.delaysetting(logical(TodayTreadmillLog.complete));
end