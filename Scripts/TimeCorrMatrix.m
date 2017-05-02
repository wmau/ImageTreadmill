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
codingCells = 'timecells';
z = true;
similarityMetric = 'corr';

switch similarityMetric 
    case 'corr'
        yLabel = 'Mean similarity (R)';
    case 'innerproduct'
        yLabel = 'Mean similarity (inner product)';
end

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
        
        load('TimeCells.mat','TodayTreadmillLog');
        nTrials{a}(s) = sum(TodayTreadmillLog.complete);
    end
end
nAllTrials = sum(cell2mat(nTrials));
    
%% Do analysis.
allR_trials = nan(5,5,nSessions);               %5x5xS matrix of across-trial correlations, binned by trial blocks.
allR_days = nan(5,5,nAnimals);                  %5x5xA matrix of across-trial correlations, binned by day.
session=1;
trial=1;
dayR = cell(nAnimals,1);
sessionTracker = nan(nSessions,1);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions = length(ssns);
    
    [R,lapNum,sessionNum] = PVTrialCorr2(fulldataset(ssns),'codingcells',codingCells,...
        'z',z,'similarityMetric',similarityMetric);
    R_allCells = nanmean(R,3);
    
    allR_trials(:,:,session:session+length(ssns)-1) = avgCorrMatrixOverAllDays(R_allCells,lapNum,...
        sessionNum,'binByNTrials');
    session = session+length(ssns);
    
    allR_days(:,:,a) = avgCorrMatrixOverAllDays(R_allCells,lapNum,...
        sessionNum,'binByDay');
    [~,dayR{a}] = binCoeffs(R_allCells,'processingMode','binByDay','sessionNum',sessionNum); 
end

A = cell(nSessions);
for a=1:nAnimals
    for s1=1:nSessions
        for s2=1:nSessions
            A{s1,s2} = [A{s1,s2}; dayR{a}{s1,s2}];
        end
    end
end
[m,sem,diags] = collapseByLag(allR_days);
            
MATRIX = nanmean(allR_days,3);
trialMATRIX = nanmean(allR_trials,3);

figure('Position',[520   125   640   675]);
subplot(2,2,1);
imagesc(MATRIX);
axis equal; axis tight; colormap hot;
set(gca,'ticklength',[0 0],'linewidth',4,'fontsize',12);
xlabel('Day','fontsize',15);
ylabel('Day','fontsize',15);
title('Bins = 1 day');

subplot(2,2,2);
errorbar(0:4,m,sem,'linewidth',4,'capsize',0);
set(gca,'tickdir','out','linewidth',4,'fontsize',12);
xlabel('Lag','fontsize',15);
ylabel(yLabel,'fontsize',15); 
xlim([-.5 4.5]);

subplot(2,2,3);
imagesc(trialMATRIX);
axis equal; axis tight; colormap hot;
set(gca,'ticklength',[0 0],'linewidth',4,'fontsize',12);
xlabel('Trial block','fontsize',15);
ylabel('Trial block','fontsize',15);
title('Bins = multiple trials');

[m,sem,diags] = collapseByLag(allR_trials);

subplot(2,2,4);
errorbar(0:4,m,sem,'linewidth',4,'capsize',0);
set(gca,'tickdir','out','linewidth',4,'fontsize',12);
xlabel('Lag','fontsize',15);
ylabel(yLabel,'fontsize',15); 
xlim([-.5 4.5]);

X = cell2mat(diags');
lags = [];
for i=1:length(diags)
    lags = [lags (i-1).*ones(1,length(diags{i}))];
end
[~,~,stats] = anovan(X,{lags},'display','off');
figure;
c = multcompare(stats);