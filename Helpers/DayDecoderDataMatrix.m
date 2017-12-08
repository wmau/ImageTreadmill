function [trainX,testX,trainY,testY,trainRuns,testRuns] = DayDecoderDataMatrix(mds,varargin)
%
%
%

%%
    %Get all the time cells and runs for each session.
    nSessions = length(mds); 
    [timeCells,randRuns,testRuns] = deal(cell(1,nSessions)); 
    for thisSession=1:nSessions
        %Get the roster of time cells.
        timeCells{thisSession} = getTimeCells(mds(thisSession)); 
        
        load(fullfile(mds(thisSession).Location,'TimeCells.mat'),...
            'ratebylap','TodayTreadmillLog');
        complete = logical(TodayTreadmillLog.complete); 
        nRuns = round(sum(complete)*.5); 
       
        %Get random runs for training data. 
        randRuns{thisSession} = randsample(1:sum(complete),nRuns)';
    end
    
    p = inputParser;
    p.addRequired('mds',@(x) isstruct(x)); 
    p.addParameter('neurons',timeCells,@(x) iscell(x));
    p.addParameter('runs',randRuns,@(x) iscell(x));
    p.addParameter('shuffle',false,@(x) islogical(x)); 
   
    p.parse(mds,varargin{:});
    neurons = p.Results.neurons;
    shuffle = p.Results.shuffle; 
    trainRuns = p.Results.runs;
    
    %Get test runs, the runs that were not included in the training set. 
    for thisSession=1:nSessions
        load(fullfile(mds(thisSession).Location,'TimeCells.mat'),...
            'TodayTreadmillLog');
        complete = logical(TodayTreadmillLog.complete); 
        
        testRuns{thisSession} = setdiff(1:sum(complete),trainRuns{thisSession});
    end    
    
%% Get cells of interest, any cell that was a time cell on at least one day.
    map = msMatchMultiSessionCells(mds,neurons); 
    map = map(all(map > 0 & ~isnan(map),2),:);                    %Must be active on all days. 

%% Build data matrices for training and test data.
    %Preallocate. 
    [trainX,testX,trainY,testY] = deal(cell(nSessions,1)); 
   
    %For each session, get the training and test data, with option to
    %shuffle. 
    for thisSession=1:nSessions
        %Input data. 
        trainX{thisSession} = reshapeRateByLap(mds(thisSession),'neurons',...
            map(:,thisSession),'runs',trainRuns{thisSession},'collapseTrials',true);
        
        testX{thisSession} = reshapeRateByLap(mds(thisSession),'neurons',...
            map(:,thisSession),'runs',testRuns{thisSession},'collapseTrials',true,...
            'shuffle',shuffle);
        
        %Response variables. 
        trainY{thisSession} = thisSession*ones(length(trainRuns{thisSession}),1);
        testY{thisSession} = thisSession*ones(length(testRuns{thisSession}),1); 
    end
    %Make into matrix for easier computation. 
    trainX = cell2mat(trainX); 
    testX = cell2mat(testX); 
    trainY = cell2mat(trainY);
    testY = cell2mat(testY); 
    
end
    