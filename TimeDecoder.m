function [Mdl,X,testX,testLaps] = TimeDecoder(md,varargin)
%
%
%

%%  
    cd(md.Location);
    load('TimeCells.mat','ratebylap','TodayTreadmillLog');
    complete = logical(TodayTreadmillLog.complete); 
    nBins = size(ratebylap,2); 
    
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('nTrainingRuns',round(sum(complete)*.5),...
        @(x) isnumeric(x)); 
    p.addParameter('neurons',getTimeCells(md),@(x) isnumeric(x)); 
    p.addParameter('shuffle',false,@(x) islogical(x)); 
    
    p.parse(md,varargin{:});
    nTrainingRuns = p.Results.nTrainingRuns;
    neurons = p.Results.neurons; 
    shuffle = p.Results.shuffle;

%% Fit the model.
    %Take a random subset of laps. This is our training data. 
    trainingLaps = randsample(sum(complete),nTrainingRuns)';
    
    %The remaining laps are our test laps. 
    testLaps = setdiff(1:sum(complete),trainingLaps);
    
    %Get predictor matrix and the test matrix. 
    X = reshapeRateByLap(md,'runs',trainingLaps,'neurons',neurons); 
    testX = reshapeRateByLap(md,'runs',testLaps,'neurons',neurons,...
        'shuffle',shuffle); 
    
    %Get response vector. 
    t = linspace(0,10,nBins)'; 
    Y = repmat(t,nTrainingRuns,1); 
    
    %Fit the model. 
    Mdl = fitcnb(X,Y,'distributionnames','mn');
    
end