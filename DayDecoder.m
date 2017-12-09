function [Mdl,trainX,testX,trainY,testY] = DayDecoder(mds,varargin)
%
%
%

%%
    nSessions = length(mds); 
    timeCells = cell(1,nSessions); 
    for thisSession=1:nSessions
        timeCells{thisSession} = getTimeCells(mds(thisSession)); 
    end
    
    p = inputParser;
    p.addRequired('mds',@(x) isstruct(x));
    p.addParameter('neurons',timeCells,@(x) iscell(x)); 
    p.addParameter('shuffle',false,@(x) islogical(x)); 
    
    p.parse(mds,varargin{:});
    neurons = p.Results.neurons;
    shuffle = p.Results.shuffle; 
    
%%
    [trainX,testX,trainY,testY,trainRuns,testRuns] = ...
        DayDecoderDataMatrix(mds,'neurons',neurons,'shuffle',shuffle); 
    
    Mdl = fitcnb(trainX,trainY,'distributionnames','mn'); 
    
end