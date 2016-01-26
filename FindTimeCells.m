function [TimeCells,ratebylap,curves,movies,T,TodayTreadmillLog] = FindTimeCells(animal,date,session,T)
%[TimeCells,ratebylap,curves,delays,x,y,time_interp] = FindTimeCells(animal,date,session,T)
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
%       ratebylap: LxBxN matrix (L=laps, B=time bins, N=neurons) depicting
%       the rate for each cell for each lap on the treadmill. 
%
%       curves: Structure array with the following fields:
%           tuning: B-length vector representing the trial mean of
%           responses.
%
%           shuffle: IxB matrix (I=number of iterations), each row is a
%           jittered tuning curve. 
%
%           sig: B-length logical vector, 1 where that bin in the tuning
%           curve is statistically larger than shuffles based on bootstrap
%           test. 
%
%           p: B-length vector, p-values for bootstrap test. 
%
%           ci: 2xB matrix, confidence intervals. First row is the upper
%           limit. 
%
%       delays: Vector depicting duration of delay per trial. 
%
%       X&Y: Aligned tracking data. 
%
%       time_interp: Aligned time vector (s).
%       
%       FT: Aligned FT.
%

%% Find time cells. 
    ChangeDirectory(animal,date,session);
    
    %Get treadmill timestamp data. 
    TodayTreadmillLog = getTodayTreadmillLog(animal,date,session);
    TodayTreadmillLog = AlignTreadmilltoTracking(TodayTreadmillLog,TodayTreadmillLog.RecordStartTime);
    
    %Get calcium imaging data. 
    load('ProcOut_minlength_3.mat','FT');
    
    %Get rate by lap matrix. 
    disp('Getting time responses for each neuron...');
    [ratebylap,x,y,time_interp,FT,TodayTreadmillLog] = getLapResponses(animal,date,session,FT,TodayTreadmillLog);  
    
    %Preallocate. 
    [nLaps,nBins,nNeurons] = size(ratebylap);
    tuningcurve = cell(nNeurons,1);
    shufflecurve = cell(nNeurons,1);
    p = cell(nNeurons,1);
    sigcurve = cell(nNeurons,1);
    ci = cell(nNeurons,1);
    
    %Include a consistency filter. Currently using arbitrary cutoff of must
    %be active for more than a quarter of the laps. 
    pLaps = 0.2;
    critLaps = round(pLaps*nLaps);
    
    %Perform permutation test on all neurons. 
    goodlaps = [];
    disp('Performing permutation test on all neurons...');
    prog = ProgressBar(nNeurons);
    for thisNeuron=1:nNeurons
        [tuningcurve{thisNeuron},shufflecurve{thisNeuron},p{thisNeuron},sigcurve{thisNeuron},ci{thisNeuron}] = ...
                TimeTuning(ratebylap(:,:,thisNeuron),TodayTreadmillLog,T);
            
        if sum(any(ratebylap(:,:,thisNeuron),2)) > critLaps
            goodlaps = [goodlaps; thisNeuron];
        end
        prog.progress;
    end
    prog.stop;
    
    %Struct array for response curves. 
    curves.tuning = tuningcurve; 
    curves.shuffle = shufflecurve;
    curves.sig = sigcurve;
    curves.p = p; 
    curves.ci = ci;
    
    %Tracking.
    movies.x = x;
    movies.y = y;
    movies.t = time_interp;
    movies.FT = FT;
    
    %Get indices of neurons that pass the test. 
    TimeCells = intersect(find(cellfun(@any,sigcurve)),goodlaps); 
    save('TimeCells.mat','TimeCells','ratebylap','curves','movies','T','TodayTreadmillLog'); 
end