%%
%Calculates the error of the trial decoder compared to shuffle. 

%%
clear;
loadMD;

fulldataset = [MD(292:303) MD(305:308)];
nSessions = length(fulldataset); 
nTrialBlocks = 8; 

animals = unique({fulldataset.Animal});
nAnimals = length(animals); 

B = 50;
[allErrors,allErrorsShuffles] = deal([]);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
    for s=1:length(ssns)
        %Get the decode error for each run. 
        [Mdl,~,testX,testLaps,trialBlockLims] = ...
            TrialDecoder(fulldataset(ssns(s)),'nTrialBlocks',nTrialBlocks); 
        [decodedTrial,postProbs] = PredictTrial(Mdl,testX,trialBlockLims,...
            testLaps,'plotit',false);
        pCorrect = TrialDecodeError(decodedTrial,testLaps,trialBlockLims); 
        
        %Take the mean across runs. 
        allErrors = [allErrors; pCorrect];
        
        %For each session, run B shuffle tests. 
        for i=1:B
            testX = reshapeRateByLap(fulldataset(ssns(s)),...
                'runs',testLaps,'collapseTrials',true);
            testX = testX(randperm(size(testX,1)),:);
            shuffleDecode = PredictTrial(Mdl,testX,trialBlockLims,...
                testLaps,'plotit',false);
            pCorrectShuffle = TrialDecodeError(shuffleDecode,testLaps,...
                trialBlockLims); 

            allErrorsShuffles = [allErrorsShuffles; pCorrectShuffle]; 
        end
    end
end

%% Plot error 
scatterBox([allErrors;allErrorsShuffles],...
    [zeros(size(allErrors));ones(size(allErrorsShuffles))],...
    'xLabels',{'Empirical','Shuffle'},'yLabel','% trial decoded correctly',...
    'boxColor',[0 0 0; 1 0 0])

p = ranksum(allErrors,allErrorsShuffles);