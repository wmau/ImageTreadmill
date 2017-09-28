function ElapsedDayTimeDecoder(trainSession,decodeSession,varargin)
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
   
    p.parse(trainSession,decodeSession,varargin{:});
    neurons = p.Results.neurons;
    nTrainingRuns = p.Results.nTrainingRuns;
    
%% 
    Mdl = TimeDecoder(trainSession,'neurons',neurons,...
        'nTrainingRuns',nTrainingRuns,'plotit',false); 
    
%%
    cd(decodeSession.Location); 
    load('TimeCells.mat','TodayTreadmillLog'); 
    decodeSessionCompleteRuns = find(TodayTreadmillLog.complete)';
    
    testX = reshapeRateByLap(decodeSession,'neurons',matchedCells(:,2),...
        'runs',decodeSessionCompleteRuns); 
    
    [decodedTime,postProbs] = PredictTime(Mdl,testX); 
    
end