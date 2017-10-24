%%
%Calculates the error of the trial decoder compared to shuffle. 

%%
clear;
loadMD;

fulldataset = [MD(292:303) MD(305:308)];
nSessions = length(fulldataset); 

animals = unique({fulldataset.Animal});
nAnimals = length(animals); 

B = 50;
[allErrors,allErrorsShuffles] = deal([]);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
    for s=1:length(ssns)
        %Get the decode error for each run. 
        [Mdl,~,testX,testLaps,trialBlockLims] = ...
            TrialDecoder(fulldataset(ssns(s))); 
        [decodedTrial,postProbs] = PredictTrial(Mdl,testX,trialBlockLims,...
            testLaps,'plotit',false);
        decodeError = TrialDecodeError(decodedTrial,testLaps,trialBlockLims); 
        
        %Take the mean across runs. 
        allErrors = [allErrors; decodeError];
        
        %For each session, run B shuffle tests. 
        for i=1:B
            testX = reshapeRateByLap(fulldataset(ssns(s)),...
                'runs',testLaps,'collapseTrials',true);
            testX = testX(randperm(size(testX,1)),:);
            shuffleDecode = PredictTrial(Mdl,testX,trialBlockLims,...
                testLaps,'plotit',false);
            decodeErrorShuffle = TrialDecodeError(shuffleDecode,testLaps,...
                trialBlockLims); 

            allErrorsShuffles = [allErrorsShuffles; decodeErrorShuffle]; 
        end
    end
end

%% Plot error 
m = [mean(allErrors) mean(allErrorsShuffles)];
sem = [standarderror(allErrors) standarderror(allErrorsShuffles)];

figure; hold on;
bar([1,1.3],m,'barwidth',.4,'facecolor','k');
errorbar(1,m(1),sem(1),'color','k');
errorbar(1.3,m(2),sem(2),'color','k');
xlim([.8 1.5]);
set(gca,'tickdir','out','xtick',[1,1.3],'xticklabel',{'Real','Shuffle'});
ylabel('Trial Block Error');
