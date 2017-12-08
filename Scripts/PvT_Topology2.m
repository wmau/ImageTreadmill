clear; close all;
clc; 

loadMD;
fulldataset = MD(292:309); 
nSessions = length(fulldataset); 
B = 500;

[TD,PD] = deal(cell(nSessions,1));
R = nan(nSessions,B); 
parfor s=1:nSessions
    disp(['Analyzing session ',num2str(s),' of ',num2str(nSessions)]);
    centroids = getNeuronCentroids(fulldataset(s)); 
    nNeurons = size(centroids,1); 
    
    TCs = getTimeCells(fulldataset(s)); 
    PCs = getPlaceCells(fulldataset(s),0.01); 
    
    TD{s} = CellPairDistances(fulldataset(s),TCs,TCs,'centroids',centroids);
    PD{s} = CellPairDistances(fulldataset(s),PCs,PCs,'centroids',centroids); 
    
    for i=1:B
        rSamp1 = randsample(nNeurons,length(TCs));
        rSamp2 = randsample(nNeurons,length(PCs)); 
    
        r = CellPairDistances(fulldataset(s),rSamp1,rSamp2,'centroids',centroids); 
        R(s,i) = mean(r(:));
    end
    
end

TD_Means = cellfun(@(x) mean(x(:)),TD);
PD_Means = cellfun(@(x) mean(x(:)),PD);
R_Means = mean(R,2); 


figure('Position',[680 535 280 440]); hold on;
for s=1:nSessions
    plot([1,2,3],[TD_Means(s), PD_Means(s), R_Means(s)],'color',[.7 .7 .7],'linewidth',2);
end
errorbar([1,2,3],[mean(TD_Means), mean(PD_Means), mean(R_Means)],...
    [std(TD_Means)/sqrt(nSessions) std(PD_Means)/sqrt(nSessions), std(R_Means)/sqrt(nSessions)],...
    'color','k','linewidth',4); 
set(gca,'tickdir','out','linewidth',4,'xtick',[1,2,3],'xticklabel',{'Time cells','Place cells','Random'},...
    'fontsize',12); 
xlim([0.5 3.5]); 
ylabel('Mean cell pair distance (microns)','fontsize',15);
