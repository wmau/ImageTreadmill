function [TimeCells,ratebylap,curves,movies,T,TodayTreadmillLog] = FindTimeCells(md,T,varargin)
%[TimeCells,ratebylap,curves,movies,T,TodayTreadmillLog] = FindTimeCells(animal,date,session,T)
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
%       MD: Session entry. 
%
%       T: Delay duration of interest. 
%
%       varargin: 'alt_input',str: Tenaspis output file. If you want to use
%       another neural data mat file.
%       
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

%% Grab inputs
    cd(md.Location);
    
    ip = inputParser;
    ip.addRequired('md',@(x) isstruct(x)); 
    ip.addRequired('T',@(x) isnumeric(x) && isscalar(x)); 
    ip.addParameter('alt_input','FinalOutput.mat',@(x) ischar(x)); 
    ip.addParameter('savename','TimeCells.mat',@(x) ischar(x)); 
    ip.parse(md,T,varargin{:});
    
    neuraldata = fullfile(md.Location,ip.Results.alt_input);
    if strcmp(ip.Results.alt_input,'ProcOut.mat'),halfwindow = 10; else halfwindow = 0; end
    savename = ip.Results.savename; 
    
%% Basic set up.
    %[~,folder] = fileparts(md.Location); 
    
    %Get treadmill timestamp data. 
    TodayTreadmillLog = getTodayTreadmillLog(md);
    TodayTreadmillLog = AlignTreadmilltoTracking(TodayTreadmillLog,TodayTreadmillLog.RecordStartTime);
    
    %Get calcium imaging data. 
    load(neuraldata,'PSAbool');
    
    %Get rate by lap matrix. 
    disp('Getting time responses for each neuron...');
    [ratebylap,x,y,time_interp,aviFrame,PSAbool,TodayTreadmillLog] = getLapResponses(md,PSAbool,TodayTreadmillLog,halfwindow);  
    
    alternation = strcmp(TodayTreadmillLog.direction,'alternation');
    
    %Preallocate. 
    [nLaps,nBins,nNeurons] = size(ratebylap);
    nComplete = sum(TodayTreadmillLog.complete);
    
    %If the mouse is on both sides of the maze, separate tuning curves for
    %left and right. 
    if alternation, nCurves = 2; 
    else            nCurves = 1; end
    
    %Preallocate. 
    tuningcurve =   cell(nNeurons,nCurves);
    shufflecurve =  cell(nNeurons,nCurves);
    p =             cell(nNeurons,nCurves);
    sigcurve =      cell(nNeurons,nCurves);
    ci =            cell(nNeurons,nCurves);
    
    %Include a consistency filter. Currently using arbitrary cutoff of must
    %be active for more than a quarter of the laps. 
    pLaps = 0.25;
    if alternation
        pLaps = pLaps*2;
        critLaps = [round(pLaps*sum(TodayTreadmillLog.choice==1 & TodayTreadmillLog.complete)),...
                    round(pLaps*sum(TodayTreadmillLog.choice==2 & TodayTreadmillLog.complete))];
    else
        critLaps = round(pLaps*nComplete);
    end
    
    %Perform permutation test on all neurons. 
    goodlaps = [];
    disp('Performing permutation test on all neurons...');
    prog = ProgressBar(nNeurons);
    
    %Alternation. 
    if alternation
        for thisNeuron=1:nNeurons
            for turn=1:2
                good = TodayTreadmillLog.choice == turn;        %Indices for treadmill runs for left and right turns. 
                [tuningcurve{thisNeuron,turn},shufflecurve{thisNeuron,turn},...
                    p{thisNeuron,turn},sigcurve{thisNeuron,turn},ci{thisNeuron,turn}]...
                    = TimeTuning(ratebylap(good,:,thisNeuron),...
                    TodayTreadmillLog.delaysetting(good),TodayTreadmillLog.complete(good),T);
                
                if sum(any(ratebylap(good,:,thisNeuron),2)) > critLaps(turn)
                    goodlaps = [goodlaps; thisNeuron];
                end
            end
        prog.progress;
        end
    else
        for thisNeuron=1:nNeurons
            [tuningcurve{thisNeuron},shufflecurve{thisNeuron},p{thisNeuron},sigcurve{thisNeuron},ci{thisNeuron}] = ...
                    TimeTuning(ratebylap(:,:,thisNeuron),TodayTreadmillLog.delaysetting,TodayTreadmillLog.complete,T);

            if sum(any(ratebylap(:,:,thisNeuron),2)) > critLaps
                goodlaps = [goodlaps; thisNeuron];
            end
            prog.progress;
        end
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
    movies.aviFrame = aviFrame;
    movies.FT = PSAbool;
    
    %Get indices of neurons that pass the test. 
    TimeCells = intersect(find(any(cellfun(@any,sigcurve),2)),goodlaps);   %Lap criterion.
    nConsecPeaks = nan(length(TimeCells),size(sigcurve,2));
    for t = 1:size(sigcurve,2)
        n = 1;
        for tc = TimeCells'
            d = diff([0 sigcurve{tc,t} 0]);
            try
                nConsecPeaks(n,t) = max(find(d<0) - find(d>0));
            catch
                nConsecPeaks(n,t) = 0;
            end

            n=n+1; 
        end
    end
    TimeCells(all(nConsecPeaks < 2, 2)) = [];
    save(savename,'TimeCells','ratebylap','curves','movies','T','TodayTreadmillLog','alternation'); 
end