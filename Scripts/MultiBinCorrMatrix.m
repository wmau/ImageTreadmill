%% Set up.
clear;

loadMD;
%MD(292:295) = G45.
%MD(296:299) = G48.
%MD(300:304) = Bellatrix.
%MD(305:309) = Polaris.
fulldataset = MD(292:309);   

%Parameters to change.
timeCellsOnly = true;
useIntervals = true;

%Number of animals and sessions.
animals = unique({fulldataset.Animal});
nAnimals = length(animals);
nSessions = length(fulldataset); 

c = 1; session = 1;
yLims = [];
figure;
for nBins=3:8
    allR_trials = nan(nBins,nBins,nSessions);
    
    for a=1:nAnimals        
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
        nSessions = length(ssns);
        
        for s=1:nSessions
            [R,p,lapNum,sessionNum] = PVTrialCorr(fulldataset(ssns),'timeCellsOnly',timeCellsOnly,...
                'useIntervals',useIntervals,'plotit',false);
            
            allR_trials(:,:,session:session+length(ssns)-1) = avgCorrMatrixOverAllDays(R,lapNum,...
                sessionNum,'binByNTrials','nBins',nBins);
            session = session+length(ssns);
        end
    end
    
    [TrialM,TrialSEM] = collapseByLag(allR_trials);
    
    if ismember(nBins,[3,4,5])
        subplot(4,3,c);
            imagesc(nanmean(allR_trials,3));
            axis equal; axis tight; colormap hot;
            set(gca,'tickdir','out');
        subplot(4,3,c+3);
            errorbar(0:nBins-1,TrialM,TrialSEM,'linewidth',4);
            set(gca,'tickdir','out');
    elseif ismember(nBins,[6,7,8])
        subplot(4,3,c+3)
            imagesc(nanmean(allR_trials,3));
            axis equal; axis tight; colormap hot;
            set(gca,'tickdir','out');
        subplot(4,3,c+6);
            errorbar(0:nBins-1,TrialM,wTrialSEM,'linewidth',4);
            set(gca,'tickdir','out');   
    end
    
    c = c+1;
end
            