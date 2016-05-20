function [TimeCells,ratebylap,curves,movies,T,TodayTreadmillLog] = FindTimeCells(MD,T,varargin)
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

%% Basic set up.
    animal = MD.Animal;
    date = MD.Date;
    session = MD.Session; 
    
    dirstr = ChangeDirectory(animal,date,session);
    [~,folder] = fileparts(dirstr); 
    
    %Get treadmill timestamp data. 
    TodayTreadmillLog = getTodayTreadmillLog(animal,date,session);
    TodayTreadmillLog = AlignTreadmilltoTracking(TodayTreadmillLog,TodayTreadmillLog.RecordStartTime);
    
    %Name of neural data file. Changed post-Tenaspis 2. 
    neuraldata = 'ProcOut.mat';         %Default: Tenaspis 1. 
    savename = 'TimeCells.mat';         %Default.
    halfwindow = 10;                    %Default: Tenaspis 1. 
    if ~isempty(varargin)
        if any(strcmp('alt_input',varargin))
            filename = varargin{find(strcmp('alt_input',varargin))+1};
            neuraldata = fullfile(pwd,varargin{find(strcmp('alt_input',varargin))+1}); 
            
            if strcmp(filename,'T2output.mat')
                halfwindow = 0;
            end
        end
        
        if any(strcmp('savename',varargin))
            savename = varargin{find(strcmp('savename',varargin))+1}; 
        end
    end
    
    %Get calcium imaging data. 
    load(neuraldata,'FT');
    
    %Get rate by lap matrix. 
    disp('Getting time responses for each neuron...');
    [ratebylap,x,y,aviFrame,FT,TodayTreadmillLog] = getLapResponses(MD,FT,TodayTreadmillLog,halfwindow);  
    
    alternation = strcmp(TodayTreadmillLog.direction,'alternation');
    blocked = ~isempty(strfind(folder,'blocked')); 
    
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
    movies.t = aviFrame;
    movies.FT = FT;
    
    %Get indices of neurons that pass the test. 
    TimeCells = intersect(find(cellfun(@any,sigcurve)),goodlaps); 
    temp = cellfun(@(a) max(diff( [0 (find( ~ (a > 0) ) )+ 1] ))-1,sigcurve,'unif',0);
    longpeak = find(cellfun(@(a) any(a>1),temp));
    TimeCells = intersect(TimeCells,longpeak);
    save(savename,'TimeCells','ratebylap','curves','movies','T','TodayTreadmillLog','alternation'); 
end