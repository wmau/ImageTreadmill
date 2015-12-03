function ratebylap = getLapResponses(animal,date,sessionNum,FT,TodayTreadmillLog)
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

%% Preliminary stuff. 
    ChangeDirectory(animal,date,sessionNum);
    
    %Align FT. 
    [~,~,~,FT,~,~,~,time_interp] = AlignImagingToTracking(0.15,FT);
    
    %Get treadmill run epochs. 
    inds = getTreadmillEpochs(TodayTreadmillLog,time_interp);
    nFramesBetween = diff(inds,1,2);        %Number of frames epochs span. 
    
%% Bin time responses. 
    %Initialize.  
    [nNeurons,nFrames] = size(FT); 
    tResolution = 0.5;                                      %seconds
    nBins = TodayTreadmillLog.delaysetting/tResolution;     %vector specifying number of bins per lap
    ratebylap = nan(TodayTreadmillLog.numRuns,max(nBins),nNeurons); 
    bin = nan(TodayTreadmillLog.numRuns,max(nBins),nNeurons); 
    
    p = ProgressBar(nNeurons);
    for thisNeuron=1:nNeurons
        for thisLap=1:TodayTreadmillLog.numRuns
            tStart = 0;
            tEnd = TodayTreadmillLog.stopts(thisLap) - TodayTreadmillLog.startts(thisLap); 
            %Time vector. 
            t = linspace(tStart,tEnd,nFramesBetween(thisLap)+1); 
           
            %Times where there was a spike. 
            tspk = t(FT(thisNeuron,inds(thisLap,1):inds(thisLap,2))==1); 
            
            %Edges for histogram.
            edges = linspace(0,TodayTreadmillLog.delaysetting(thisLap),nBins(thisLap));
            
            %Make histogram.
            [ratebylap(thisLap,1:nBins(thisLap),thisNeuron),...
                bin(thisLap,1:nBins(thisLap),thisNeuron)] = hist(tspk,edges); 
            
        end
        p.progress;
    end
    p.stop;

    
end