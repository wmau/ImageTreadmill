function [decodedTime,postProbs,Mdl] = ...
    ElapsedDayTimeDecoder(trainSession,decodeSession,varargin)
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
    p.addParameter('nTrainingRuns',nCompleteRuns,@(x) isnumeric(x)); 
    p.addParameter('shuffle',false,@(x) islogical(x)); 
    p.addParameter('plotit',true,@(x) islogical(x)); 
    p.addParameter('Mdl',[]);
   
    p.parse(trainSession,decodeSession,varargin{:});
    neurons = p.Results.neurons;
    nTrainingRuns = p.Results.nTrainingRuns;
    shuffle = p.Results.shuffle;
    plotit = p.Results.plotit;
    Mdl = p.Results.Mdl;
    
%% 
    if isempty(Mdl)
        Mdl = TimeDecoder(trainSession,'neurons',neurons,...
            'nTrainingRuns',nTrainingRuns); 
    end
    
%%
    cd(decodeSession.Location); 
    load('TimeCells.mat','TodayTreadmillLog'); 
    decodeSessionCompleteRuns = find(TodayTreadmillLog.complete)';
    
    testX = reshapeRateByLap(decodeSession,'neurons',matchedCells(:,2),...
        'runs',decodeSessionCompleteRuns,'shuffle',shuffle); 
    
    [decodedTime,postProbs] = PredictTime(Mdl,testX,'plotit',plotit); 
    
end