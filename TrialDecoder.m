function [Mdl,X,testX,testLaps,trialBlockLims] = TrialDecoder(md,varargin)

%%  
    cd(md.Location);
    load('TimeCells.mat','ratebylap','TodayTreadmillLog');
    complete = logical(TodayTreadmillLog.complete); 
    nBins = size(ratebylap,2); 
    
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('neurons',getTimeCells(md),@(x) isnumeric(x)); 
    p.addParameter('shuffle',false,@(x) islogical(x)); 
    p.addParameter('nTrialBlocks',6,@(x) isnumeric(x)); 
    
    p.parse(md,varargin{:});
    neurons = p.Results.neurons; 
    shuffle = p.Results.shuffle;
    nTrialBlocks = p.Results.nTrialBlocks;

%% Fit the model.
    %Get the bin limits for trial blocks. 
    trialBlockLims = linspace(0,sum(complete),nTrialBlocks+1); 
    trialBlockLims = floor(trialBlockLims);
    blockSizes = diff(trialBlockLims);
    [Y,trainingLaps] = deal([]);
    
    for b = 1:nTrialBlocks
        block = trialBlockLims(b)+1:trialBlockLims(b+1); 
        sampleTrainingLaps = randsample(block,floor(blockSizes(b)/2)); 
        trainingLaps = [trainingLaps, sampleTrainingLaps];
        
        Y = [Y; b*ones(floor(blockSizes(b)/2),1)];
    end

    %The remaining laps are our test laps. 
    testLaps = setdiff(1:sum(complete),trainingLaps);
    
    %Get predictor matrix and the test matrix. 
    X = reshapeRateByLap(md,'runs',trainingLaps,'neurons',neurons,...
        'collapseTrials',true); 
    testX = reshapeRateByLap(md,'runs',testLaps,'neurons',neurons,...
        'collapseTrials',true,'shuffle',shuffle); 
     
    %Fit the model. 
    Mdl = fitcnb(X,Y,'distributionnames','mn');
    
end