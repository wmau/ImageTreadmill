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
            'nTrainingRuns',nTrainingRuns,'shuffle',shuffle); 
    end
    
%%
    cd(decodeSession.Location); 
    load('TimeCells.mat','TodayTreadmillLog'); 
    decodeSessionCompleteRuns = 1:sum(TodayTreadmillLog.complete)';
    
    day2Cells = matchedCells(ismember(matchedCells(:,1),neurons),2);
    testX = reshapeRateByLap(decodeSession,'neurons',day2Cells,...
        'runs',decodeSessionCompleteRuns); 
    
    [decodedTime,postProbs] = PredictTime(Mdl,testX,'plotit',plotit); 
    
end