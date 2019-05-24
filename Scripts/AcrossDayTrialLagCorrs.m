clear; loadMD;
fulldataset = [MD(292:303) MD(305:308)]; 

animals = unique({fulldataset.Animal});
nAnimals = length(animals);

figure('Position',[351    61   984   879]);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    thisData = fulldataset(ssns);
    
    [R,lapNum,sessionNum] = PVTrialCorr2(thisData);
    nDays = length(thisData); 

    maxLag = max(lapNum-1);
    binnedRs = cell(maxLag+1,nDays-1);

    dayCombs = combnk(1:nDays,2); 
    nDayCombs = size(dayCombs,1); 
    for thisDayComb = 1:nDayCombs
        referenceDay = dayCombs(thisDayComb,1); 
        comparisonDay = dayCombs(thisDayComb,2); 

        nTrials = max(lapNum(sessionNum==referenceDay)); 

        trialCombs = combnk(1:nTrials,2); 
        trialCombs = [trialCombs; repmat((1:nTrials)',1,2)];
        %trialCombs = [ones(nTrials,1), (1:nTrials)'];

        nTrialCombs = size(trialCombs,1); 

        dayDifference = comparisonDay-referenceDay;

        for thisTrialComb = 1:nTrialCombs
            referenceTrial = trialCombs(thisTrialComb,1); 
            comparisonTrial = trialCombs(thisTrialComb,2); 

            row = (lapNum==referenceTrial) & (sessionNum==referenceDay);
            col = (lapNum==comparisonTrial) & (sessionNum==comparisonDay);

            trialDifference = comparisonTrial-referenceTrial;        
            insertMe = squeeze(R(row,col,:));

            binnedRs{trialDifference+1,dayDifference} = ...
                [binnedRs{trialDifference+1,dayDifference};...
                insertMe];
        end
    end


    m = cellfun(@nanmean,binnedRs,'unif',0);
    sem = cellfun(@standarderror,binnedRs,'unif',0); 

    good = all(~cellfun(@(x) isempty(x),binnedRs),2);
    m = m(good,:);
    sem = sem(good,:);

    m = cell2mat(m);
    sem = cell2mat(sem); 

    flatM = m(:);
    flatSEM = sem(:);

    subplot(2,2,a);
    errorbar(flatM,flatSEM,'linewidth',2);
    yLims = get(gca,'ylim');
    line([size(m,1) size(m,1)],yLims,'color','r','linewidth',2);
    line([size(m,1)*2 size(m,1)*2],yLims,'color','r','linewidth',2);
    axis tight;
    xlabel('Lag');
    ylabel('R');
    title(animals{a});
end