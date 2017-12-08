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
z = false;
slidingWindow = true;
similarityMetric = 'corr';

switch similarityMetric
    case 'corr'
        yAxisLabel = 'Corr Coeff';
    case 'innerproduct'
        yAxisLabel = 'Inner Product';
end

%Number of animals and sessions.
animals = unique({fulldataset.Animal});
nAnimals = length(animals);
nSessions = length(fulldataset); 


c = 1; session = 1;
nBins = 5;
figure('Position',[-1215 -60 455 975]);
binSizeVector = 1:6;
[cMin,cMax] = deal(zeros(length(binSizeVector),1));
yLims = [];
for binSize=binSizeVector
    matrix = cell(nBins,nBins,nSessions);
    
    for a=1:nAnimals        
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
        nSessions = length(ssns);
        
        [R,lapNum,sessionNum] = PVTrialCorr2(fulldataset(ssns),...
            'codingCells',codingCells,'similarityMetric',similarityMetric,...
            'z',z);
        R_allCells = nanmean(R,3);
        
        for s=1:nSessions
            dayR = R_allCells(sessionNum==s,sessionNum==s); 
            matrix(:,:,ssns(s)) = weirdassBinningMethod(dayR,binSize,nBins,slidingWindow);
            
            %matrix(:,:,ssns(s)) = cellfun(@nanmean,temp);
        end
    end
    
    [TrialM,TrialSEM] = collapseByLag(matrix);
    %MATRIX = nanmean(matrix,3);
    MATRIX = nanmean(cellfun(@nanmean,matrix),3);
    cMin(binSize) = min(MATRIX(:));
    cMax(binSize) = max(MATRIX(:));
    
    if ismember(binSize,[1,2,3])
        subplot(4,3,c);
            imagesc(MATRIX);
            axis equal; axis tight; colormap hot;           
            set(gca,'tickdir','out');
            title(['Bin size = ',num2str(binSize)]);
            if binSize==1
                xlabel('Trial block');
                ylabel('Trial block');
            end
        subplot(4,3,c+3);
            errorbar(0:nBins-1,TrialM,TrialSEM,'linewidth',4);
            set(gca,'tickdir','out');   axis tight;
            yLims = [yLims, get(gca,'ylim')];
            if binSize==1
                ylabel(yAxisLabel);
            end
    elseif ismember(binSize,[4,5,6])
        subplot(4,3,c+3)
            imagesc(MATRIX);
            axis equal; axis tight; colormap hot;
            set(gca,'tickdir','out');
            title(['Bin size = ',num2str(binSize)]);
        subplot(4,3,c+6);
            errorbar(0:nBins-1,TrialM,TrialSEM,'linewidth',4);
            set(gca,'tickdir','out');   axis tight;
            yLims = [yLims, get(gca,'ylim')];
            if binSize==5
                xlabel('Lags');
            end        
    end
    
    c = c+1;
end

CMAX = max(cMax);
CMIN = min(cMin);

YMAX = max(yLims);
YMIN = min(yLims);

for i=binSizeVector
    if ismember(i,[1,2,3])
        c = 0; d = 3;
    elseif ismember(i,[4,5,6])
        c = 3; d = 6;
    end 
    subplot(4,3,i+c);
    caxis([CMIN CMAX]);
    
    subplot(4,3,i+d);
    ylim([YMIN YMAX]);
end
            