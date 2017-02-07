clear
loadMD;
fulldataset = MD(292:309); 
animals = unique({fulldataset.Animal});
nAnimals = length(animals);
colors = parula(nAnimals);

statType = 'ti';
switch statType
    case 'ti', c = [0 .5 .5];
    case 'si', c = [.58 .44 .86];
end
PCcrit = .01;

S = []; US = [];
[SM,USM] = deal(cell(nAnimals,1));
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions = length(ssns)-1;
    [SM{a},USM{a}] = deal(cell(nSessions,1));
    
    for s=1:nSessions
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
        SM{a}{s} = mean(temp(stable,:));
        US = [US; temp(unstable,:)];
        USM{a}{s} = mean(temp(unstable,:));
    end
    
end

%dpvalue = signrank(s(:,1),s(:,2));
% [sR,sp] = corr(S(:,1),S(:,2),'type','spearman');
% [usR,usp] = corr(US(:,1),US(:,2),'type','spearman');

figure('Position',[260 270 560 420]); hold on;
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

SM = cell2mat(cellfun(@cell2mat, SM,'unif',0));
USM = cell2mat(cellfun(@cell2mat, USM,'unif',0));
nSessions = length(SM);
SSEM = std(SM)./sqrt(nSessions);
USSEM = std(USM)./sqrt(nSessions);

figure('Position',[840 230 290 470]); hold on;
    plot([1,2],[SM(:,1) SM(:,2)],'color',[c .7],'linewidth',.2);
    plot([1,2],[USM(:,1) USM(:,2)],'linestyle','--','color',[.6 .6 .6 .7],'linewidth',.2);
    errorbar([1,2],mean(SM),SSEM,'color',c,'linewidth',5);
    errorbar([1,2],mean(USM),USSEM,'color',[.6 .6 .6],'linewidth',5);

ylabel('Mean Mutual Information [bits]');
set(gca,'xticklabel',{'Day 1','Day 2'});
set(gca,'xtick',[1:2],'tickdir','out');
xlim([0.5,2.5]);
legend({'Stable','Unstable'},'Location','southwest');

stability = [ones(nSessions*2,1); zeros(nSessions*2,1)];
session = repmat([ones(nSessions,1); 2*ones(nSessions,1)],[2,1]);
grps = {stability,session};
X = [SM(:);USM(:)];
[p,tbl,stats,terms] = anovan(X,grps,'model','full','varnames',{'Stability','Session'},...
    'display','off');
comps = multcompare(stats,'dimension',[1,2],'display','off');
disp(['Main effect of stability F = ',num2str(tbl{2,6}),', P = ',num2str(tbl{2,7})]);
disp(['Main effect of session F = ',num2str(tbl{3,6}),', P = ',num2str(tbl{3,7})]);
disp(['Interaction F = ',num2str(tbl{4,6}),', P = ',num2str(tbl{4,7})]);
disp(['Difference between Day 1 and Day 2 for unstable cells P = ',num2str(comps(2,6))]);
disp(['Difference between Day 1 and Day 2 for stable cells P = ',num2str(comps(5,6))]);
%%
% nS = size(S,1); nUS = size(US,1); 
% stability = [ones(nS*2,1); zeros(nUS*2,1)];
% session = [ones(nS,1); 2*ones(nS,1); ones(nUS,1); 2*ones(nUS,1)];
% X = [S(:); US(:)];
% subjects = [repmat([1:nS]',[2,1]); repmat([1:nUS]',[2,1])];
% grps = {stability,session};
% [p,tbl,stats,terms] = anovan(X,grps,'model','full');

% figure; hold on;
% sM = mean(S); 
% usM = mean(US); 
% sSEM = std(S)./sqrt(nS);
% usSEM = std(US)./sqrt(nUS);
% figure; hold on;
% errorbar([1,2],sM,sSEM);
% errorbar([1.1,2.1],usM,usSEM);