%% Description
% Performs pairwise correlations of population vectors across scales. 

%% Set up.
clear;

loadMD;
%MD(292:295) = G45.
%MD(296:299) = G48.
%MD(300:304) = Bellatrix.
%MD(305:309) = Polaris.
fulldataset = MD(292:309);   

%Parameters to change.
placeCellsOnly = true;

%Number of animals and sessions.
animals = unique({fulldataset.Animal});
nAnimals = length(animals);
nSessions = length(fulldataset); 

%Get number of trials. 
nTrials = cell(nAnimals,1);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions = length(ssns);
    
    nTrials{a} = zeros(nSessions,1);
    for s=1:nSessions
        cd(fulldataset(ssns(s)).Location); 
        
        load('Alternation.mat','Alt');
        nTrials{a}(s) = size(Alt.summary,1);
    end
end
nAllTrials = sum(cell2mat(nTrials));
    
%% Do analysis.
allR_trials = nan(5,5,nSessions);               %5x5xS matrix of across-trial correlations, binned by trial blocks.
allR_days = nan(5,5,nAnimals);                  %5x5xA matrix of across-trial correlations, binned by day.
session=1;
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions = length(ssns);
  
    [R,lapNum,sessionNum] = PVTrialCorr_place(fulldataset(ssns),'placeCellsOnly',placeCellsOnly);
    
    allR_trials(:,:,session:session+length(ssns)-1) = avgCorrMatrixOverAllDays(R,lapNum,...
        sessionNum,'binByNTrials');
    session = session+length(ssns);
    
    allR_days(:,:,a) = avgCorrMatrixOverAllDays(R,lapNum,sessionNum,'binByDay');
end

close all;
figure;
subplot(2,2,1);
imagesc(nanmean(allR_trials,3));
axis equal; axis tight; colormap hot;
set(gca,'tickdir','out');
xlabel('Trial Blocks');
ylabel('Trial Blocks');
title('Bins = multiple trials');

subplot(2,2,2);
imagesc(nanmean(allR_days,3));
axis equal; axis tight; colormap hot;
set(gca,'tickdir','out');
xlabel('Day');
ylabel('Day');
title('Bins = 1 day');

subplot(2,2,3);
[TrialM,TrialSEM] = collapseByLag(allR_trials);
errorbar(0:4,TrialM,TrialSEM,'linewidth',4);
set(gca,'tickdir','out');
xlabel('Lag');

subplot(2,2,4);
[dayM,daySEM] = collapseByLag(allR_days);
errorbar(0:4,dayM,daySEM,'linewidth',4);
set(gca,'tickdir','out');
xlabel('Lag');
