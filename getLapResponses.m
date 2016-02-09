function [ratebylap,x,y,aviFrame,FT,TodayTreadmillLog] = getLapResponses(animal,date,sessionNum,FT,TodayTreadmillLog)
%[ratebylap,delays,x,y,time_interp,FT] = getLapResponses(animal,date,sessionNum,FT,TodayTreadmillLog)
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
%       X&Y: Aligned position data. 
%       
%       time_interp: Interpolated and aligned timestamps. 
%
%       FT: From the input, but aligned. 
%
%       TodayTreadmillLog: Containing corrections in .complete for where
%       mouse is on treadmill. 
%

%% Preliminary stuff. 
    dirstr = ChangeDirectory(animal,date,sessionNum);
    [~,folder] = fileparts(dirstr); 
    alternation = strcmp(TodayTreadmillLog.direction,'alternation'); 
    blocked = ~isempty(strfind(folder,'blocked')); 
    
    %Align FT. 
    [x,y,~,FT,~,~,aviFrame] = AlignImagingToTracking(0.15,FT);
    
    %Get treadmill run epochs. 
    inds = getTreadmillEpochs(TodayTreadmillLog,aviFrame);
    nFramesBetween = diff(inds,1,2);                        %Number of frames epochs span. 
    
    %Make sure the mouse is on the treadmill during these epochs. Find the
    %section that the mouse is on. 
    if alternation
        Alt = postrials_treadmill(x,y,0); 
        if blocked %In blocked alternation, choice is always correct.
            Alt.alt = ones(1,length(Alt.frames)); 
        end
        sect = Alt.section; 
    else
        bounds = sections_treadmill(x,y,TodayTreadmillLog.direction);
        sect = getsection_treadmill(x,y,bounds);
    end
    
%% Bin time responses. 
    %Initialize.  
    [nNeurons,~] = size(FT); 
    tResolution = 0.25;                                     %seconds
    nBins = TodayTreadmillLog.delaysetting/tResolution;     %Vector specifying number of bins per lap
    nRuns = TodayTreadmillLog.numRuns;
    
    %Preallocate.
    TodayTreadmillLog.choice = zeros(nRuns,1); 
    
    %Preallocate.
    ratebylap = nan(nRuns,max(nBins),nNeurons);         
    propOnTM = nan(nRuns,1);
    onTM = nan(nRuns,1);
    
    for lapnum=1:nRuns       
        %Proportion of time spent on treadmill. 
        propOnTM(lapnum) = sum(sect(inds(lapnum,1):inds(lapnum,2))==2)/(nFramesBetween(lapnum)+1);
        
        %Logical. Currently using threshold that half the time must be
        %spent on the treadmill.
        onTM(lapnum) = propOnTM(lapnum) > 0.5;
        
        if ~onTM(lapnum)
            TodayTreadmillLog.complete(lapnum) = 0;
        end
        
        if TodayTreadmillLog.complete(lapnum) && alternation
            %Left or right. 
            TodayTreadmillLog.choice(lapnum) = mode(Alt.choice(inds(lapnum,1):inds(lapnum,2))); 
        end
    end
    
    p = ProgressBar(nNeurons);
    for thisNeuron=1:nNeurons
        for lapnum=1:nRuns           
            if TodayTreadmillLog.complete(lapnum)
                tStart = 0;
                tEnd = TodayTreadmillLog.stopts(lapnum) - TodayTreadmillLog.startts(lapnum); 
                %Time vector. 
                t = linspace(tStart,tEnd,nFramesBetween(lapnum)+1); 
                
                %Times where there was a spike. 
                tspk = t(FT(thisNeuron,inds(lapnum,1):inds(lapnum,2))==1);

                %Edges for histogram.
                edges = linspace(0,TodayTreadmillLog.delaysetting(lapnum),nBins(lapnum));

                %Make histogram.
                [ratebylap(lapnum,1:nBins(lapnum),thisNeuron),~] = hist(tspk,edges);  
            end
        end
        p.progress;
       
    end
    p.stop;

end