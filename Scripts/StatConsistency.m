clear
loadMD;
fulldataset = MD(292:309); 
animals = unique({fulldataset.Animal});
nAnimals = length(animals);
colors = parula(nAnimals);

statType = 'si';
switch statType
    case 'ti', c = [0 .5 .5];
    case 'si', c = [.58 .44 .86];
end
PCcrit = .01;

S = []; US = [];
for a=1:nAnimals
    N = 0;
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
    for s=1:length(ssns)-1
        cd(fulldataset(ssns(s)).Location);
        load('FinalOutput.mat','NumNeurons');
        
        switch statType
            case 'ti',neurons = getTimeCells(fulldataset(ssns(s))); 
            case 'si',neurons = getPlaceCells(fulldataset(ssns(s)),PCcrit);
        end
        
        neurons = EliminateUncertainMatches([fulldataset(ssns(s)) fulldataset(ssns(s+1))],...
            neurons);
        
        switch statType
            case 'si', corrs = CorrPlaceFields(fulldataset(ssns(s)),fulldataset(ssns(s+1)),neurons);
            case 'ti', corrs = CorrTrdmllTrace(fulldataset(ssns(s)),fulldataset(ssns(s+1)),neurons);
        end
        
        stblcrit = .01/length(neurons);
        stable = intersect(find(corrs(:,2) < stblcrit),neurons);
        unstable = intersect(find(corrs(:,2) > stblcrit | isnan(corrs(:,2))),neurons);
        
        temp = msStats(fulldataset(ssns(s):ssns(s+1)),statType,1:NumNeurons);
        S = [S; temp(stable,:)];
        US = [US; temp(unstable,:)];
        
        N = N + size(temp,1); 
    end
    
end

%dpvalue = signrank(s(:,1),s(:,2));
% [sR,sp] = corr(S(:,1),S(:,2),'type','spearman');
% [usR,usp] = corr(US(:,1),US(:,2),'type','spearman');

figure; hold on;
scat = scatter(S(:,1),S(:,2),20,c,'filled');
alpha(scat,.5);
scat = scatter(US(:,1),US(:,2),20,'k','filled');
alpha(scat,.5);
l = lsline;
[l(1:2).LineStyle] = deal('--'); 
[l(1:2).LineWidth] = deal(3);
l(1).Color = [.7 .7 .7];
l(2).Color = c;
set(gca,'tickdir','out');
axis tight;
yLims = get(gca,'ylim'); xLims = get(gca,'xlim');
line(xLims,yLims,'color','k','linewidth',2);
xlabel('Day 1'); ylabel('Day 2');
title('Mutual Spatial Information [bits]');
%text(xLims(1),yLims(2),['Stable P = ',num2str(sR), ' Unstable P = ',num2str(usR)]);

sM = mean(S);
usM = mean(US); 
figure; hold on;