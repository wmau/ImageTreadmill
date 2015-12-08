function TimeCells = FindTimeCells(sessionStruct,animal,date,session,T)
%TimeCells = FindTimeCells(sessionStruct,animal,date,session,T)
%
%   Finds time cells using a few criteria. First, the neuron must be active
%   for at least some proportion of the laps. This proportion is hard-coded
%   in and is currently 0.25. Next, we see what happens if we temporally
%   shuffle the responses. For each iteration, we circularly permute each
%   lap's responses and take the mean across laps to produce a surrogate
%   tuning curve. Then we compare the empirical tuning curve and look for
%   regions where the empirical tuning curve exceeds the shuffled data. A
%   time cell is declared if it passes both these tests. 
%
%   INPUTS
%       sessionStruct: MD entry. 
%
%       animal: Mouse (e.g., GCamp6f_45_treadmill).
%
%       date: Recording date (e.g., 11_19_2015).
%
%       session: Session number. 
%
%       T: Delay duration of interest. 
%
%   OUTPUT
%       TimeCells: Index referencing FT of neurons that pass the tests. 
%

%%
    cd(sessionStruct.Location);
    
    TodayTreadmillLog = getTodayTreadmillLog(animal,date,session);
    TodayTreadmillLog = AlignTreadmilltoTracking(TodayTreadmillLog,TodayTreadmillLog.RecordStartTime);
    
    load('ProcOut.mat','FT');
    
    %Get rate by lap matrix. 
    disp('Getting time responses for each neuron...');
    [ratebylap,delays,x,y,time_interp] = getLapResponses(animal,date,2,FT,TodayTreadmillLog);  
    
    %Preallocate. 
    [nLaps,nBins,nNeurons] = size(ratebylap);
    tuningcurve = cell(nNeurons,1);
    shufflecurve = cell(nNeurons,1);
    p = cell(nNeurons,1);
    sigcurve = cell(nNeurons,1);
    
    %Include a consistency filter. Currently using arbitrary cutoff of must
    %be active for more than a quarter of the laps. 
    pLaps = 0.25;
    critLaps = pLaps*nLaps;
    
    %Perform permutation test on all neurons. 
    disp('Performing permutation test on all neurons...');
    prog = ProgressBar(nNeurons);
    for thisNeuron=1:nNeurons
        if sum(any(ratebylap(:,:,thisNeuron),2)) > critLaps
            [tuningcurve{thisNeuron},shufflecurve{thisNeuron},p{thisNeuron},sigcurve{thisNeuron}] = ...
                TimeTuning(ratebylap(:,:,thisNeuron),delays,T);
        end
        prog.progress;
    end
    prog.stop;
    
    %Get indices of neurons that pass the test. 
    TimeCells = find(cellfun(@any,sigcurve)); 

end