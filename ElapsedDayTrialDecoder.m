function [decodedTrial,postProbs,Mdl,trialBlockLims] = ...
    ElapsedDayTrialDecoder(trainSession,decodeSession,varargin)
%
%
%

%%
    cd(trainSession.Location);
    load('TimeCells.mat','TodayTreadmillLog');
    nCompleteRuns = sum(TodayTreadmillLog.complete); 
    
    trainingTCs = getTimeCells(trainSession); 
    matchedCells = msMatchCells([trainSession decodeSession],trainingTCs,true);
   
    p = inputParser;
    p.addRequired('trainSession',@(x) isstruct(x));
    p.addRequired('decodeSession',@(x) isstruct(x)); 
    p.addParameter('neurons',matchedCells(:,1),@(x) isnumeric(x)); 
    p.addParameter('pBlockforTraining',1,@(x) isnumeric(x)); 
    p.addParameter('nTrialBlocks',6,@(x) isnumeric(x)); 
    p.addParameter('shuffle',false,@(x) islogical(x)); 
    p.addParameter('plotit',true,@(x) islogical(x)); 
    p.addParameter('Mdl',[]);
   
    p.parse(trainSession,decodeSession,varargin{:});
    neurons = p.Results.neurons;
    pBlockforTraining = p.Results.pBlockforTraining;
    nTrialBlocks = p.Results.nTrialBlocks;
    shuffle = p.Results.shuffle;
    plotit = p.Results.plotit;
    Mdl = p.Results.Mdl;
    
%% 
    if isempty(Mdl)
        Mdl = TrialDecoder(trainSession,'neurons',neurons,...
            'pBlockforTraining',pBlockforTraining);
    end
    
%% 
    cd(decodeSession.Location);
    
    load('TimeCells.mat','TodayTreadmillLog'); 
    decodeSessionCompleteRuns = 1:sum(TodayTreadmillLog.complete);
    
    testX = reshapeRateByLap(decodeSession,'neurons',matchedCells(:,2),...
        'runs',decodeSessionCompleteRuns,'shuffle',shuffle,'collapseTrials',true); 
    
    trialBlockLims = linspace(0,sum(TodayTreadmillLog.complete),nTrialBlocks+1); 
    trialBlockLims = floor(trialBlockLims);
%%
    [decodedTrial,postProbs] = PredictTrial(Mdl,testX,...
        trialBlockLims,decodeSessionCompleteRuns,'plotit',plotit);
end