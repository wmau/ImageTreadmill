%% Description
% Performs pairwise correlations of population vectors across scales. 
% Does a cell-by-cell trace correlation then averages across cells and then
% across time in two ways: 
%
% (1) Dump all correlation coefficients into one day and average. 
% (2) Average across all trials. 

%% Set up.
clear;

loadMD;
%MD(292:295) = G45.
%MD(296:299) = G48.
%MD(300:304) = Bellatrix.
%MD(305:309) = Polaris.
fulldataset = [MD(292:299) MD(300:303) MD(305:308)];   

%Parameters to change.
codingCells = 'timecells';      %Options: timecells or placecells
z = true;                       
similarityMetric = 'corr';      %Options: corr or innerproduct.

switch codingCells
    case 'timecells'
        c = [0 .5 .5];
    case 'placecells'
        c = [.58 .44 .86];
end
        

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
nBins = 4;
binnedByTrials = nan(nBins,nBins,nSessions);               %5x5xS matrix of across-trial correlations, binned by trial blocks.
binnedByDays = nan(nBins,nBins,nAnimals);                  %5x5xA matrix of across-trial correlations, binned by day.
session=1;
trial=1;
dayR = cell(nAnimals,1);
R_DayBlocks = [];
sessionTracker = nan(nSessions,1);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
    %Do activity correlations then take the mean across neurons. R_allcells
    %describes the entire population.
    [allTrialRs,lapNum,sessionNum] = PVTrialCorr2(fulldataset(ssns),'codingcells',codingCells,...
        'z',z,'similarityMetric',similarityMetric);
    R_allCells = nanmean(allTrialRs,3);
    
    %Bin trials into 5 bins. 
    [binnedByTrials(:,:,session:session+length(ssns)-1),sessionRs] = ...
        avgCorrMatrixOverAllDays(R_allCells,lapNum,sessionNum,'binByNTrials','nBins',nBins);
    session = session+length(ssns);
    
    %Also collect the same-day correlations.
    R_DayBlocks = [R_DayBlocks; sessionRs]; 
    
    %Bin correlations by days. 
    binnedByDays(:,:,a) = avgCorrMatrixOverAllDays(R_allCells,lapNum,...
        sessionNum,'binByDay','nBins',nBins);
    %dayR is a cell array with one entry per animal. In dayR are 5x5 cells
    %containing the correlation values for that day by day correlation.
    [~,dayR{a}] = binCoeffs(R_allCells,'processingMode','binByDay','sessionNum',sessionNum); 
end

    matSize = min(cellfun('length',R_DayBlocks));
    
    nSessions = length(R_DayBlocks);
    allTrialRs = zeros(matSize,matSize,nSessions);
    for s=1:nSessions
        allTrialRs(:,:,s) = R_DayBlocks{s}(1:matSize,1:matSize);
    end
    
% A = cell(nSessions);
% for a=1:nAnimals
%     for s1=1:nSessions
%         for s2=1:nSessions
%             A{s1,s2} = [A{s1,s2}; dayR{a}{s1,s2}];
%         end
%     end
% end
            
daysMATRIX = nanmean(binnedByDays,3);
binTrialMATRIX = nanmean(binnedByTrials,3);
allTrialMATRIX = nanmean(allTrialRs,3);

figure('Position',[520   125   640   675]);
subplot(2,2,1);
imagesc(daysMATRIX);
axis equal; axis tight; colormap hot;
set(gca,'ticklength',[0 0],'linewidth',4,'fontsize',12);
xlabel('Day','fontsize',15);
ylabel('Day','fontsize',15);
title('Bins = 1 day');
colorbar;

[m,sem,diags] = collapseByLag(binnedByDays);
subplot(2,2,2);
errorbar(0:nBins-1,m,sem,'linewidth',4,'capsize',0,'color',c);
set(gca,'tickdir','out','linewidth',4,'fontsize',12);
xlabel('Lag','fontsize',15);
ylabel(yLabel,'fontsize',15); 
xlim([-.5 nBins-1+.5]);

subplot(2,2,3);
imagesc(allTrialMATRIX);
axis equal; axis tight; colormap hot;
set(gca,'ticklength',[0 0],'linewidth',4,'fontsize',12);
xlabel('Trial','fontsize',15);
ylabel('Trial','fontsize',15);
colorbar;

[m,sem,diags] = collapseByLag(allTrialRs);
subplot(2,2,4);
errorbar(0:matSize-1,m,sem,'linewidth',4,'capsize',0,'color',c);
set(gca,'tickdir','out','linewidth',4,'fontsize',12);
xlabel('Lag','fontsize',15);
ylabel(yLabel,'fontsize',15); 
xlim([-.5 matSize+.5]);

%% ANOVA by trials. 
% X = cell2mat(diags');
% lags = [];
% for i=1:length(diags)
%     lags = [lags (i-1).*ones(1,length(diags{i}))];
% end
%[~,~,stats] = anovan(X,{lags},'display','off');


%% ANOVA by days, all cell entries
X = [];
lags = [];
for a=1:nAnimals
    [m,sem,diags] = collapseByLag(dayR{a});
    X = [X cell2mat(diags')];
    
    for i=1:length(diags)
        lags = [lags (i-1).*ones(1,length(diags{i}))];
    end
    
end

[~,~,stats] = anovan(X,{lags},'display','off');
figure;
comp = multcompare(stats);