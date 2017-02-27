clear
loadMD;
fulldataset = MD(292:309); 
animals = unique({fulldataset.Animal});
nAnimals = length(animals);
colors = parula(nAnimals);

statType = 'ti';
switch statType
    case 'ti', col = [0 .5 .5];
    case 'si', col = [.58 .44 .86];
end
PCcrit = .01;

saveBool = true;
folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures';
scattername = fullfile(folder,['Stable ',statType,' Scatter']);
linename = fullfile(folder,['Stable ',statType,' Drop']);

if saveBool
    c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');
    
    if ~strcmp(c,'y')
        saveBool = false;
    end
end

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
        US = [US; temp(unstable,:)];
        
        if length(stable) > 4
            sTemp = temp(stable,:); 
            sTemp(sTemp(:,2)==0,:) = [];
            SM{a}{s} = mean(log10(sTemp),1);
        else, SM{a}{s} = nan(1,2);
        end
        
        if length(unstable) > 4
            uTemp = temp(unstable,:);
            uTemp(uTemp(:,2)==0,:) = [];
            USM{a}{s} = mean(log10(uTemp),1);
        else, USM{a}{s} = nan(1,2); 
        end
        
        
    end
    
end

%dpvalue = signrank(s(:,1),s(:,2));
% [sR,sp] = corr(S(:,1),S(:,2),'type','spearman');
% [usR,usp] = corr(US(:,1),US(:,2),'type','spearman');

figure('Position',[260 270 560 420]); hold on;
scat = scatter(S(:,1),S(:,2),20,col,'filled');
alpha(scat,.5);
scat = scatter(US(:,1),US(:,2),20,'k','filled');
alpha(scat,.5);
l = lsline;
[l(1:2).LineStyle] = deal('--'); 
[l(1:2).LineWidth] = deal(3);
l(1).Color = [.7 .7 .7];
l(2).Color = col;
set(gca,'tickdir','out');
axis tight;
yLims = get(gca,'ylim'); xLims = get(gca,'xlim');
line(xLims,yLims,'color','k','linewidth',2);
xlabel('Day 0'); ylabel('Day 1');
title('Mutual Info. [bits]');
%text(xLims(1),yLims(2),['Stable P = ',num2str(sR), ' Unstable P = ',num2str(usR)]);
if saveBool, print(scattername,'-dpdf'); end

SM = cell2mat(cellfun(@cell2mat, SM,'unif',0));
USM = cell2mat(cellfun(@cell2mat, USM,'unif',0));
nSessions = length(SM);
SSEM = nanstd(SM)./sqrt(nSessions);
USSEM = nanstd(USM)./sqrt(nSessions);

figure('Position',[840 230 290 470]); hold on;
    plot([1,2],[SM(:,1) SM(:,2)],'color',[col .7],'linewidth',.2);
    plot([1,2],[USM(:,1) USM(:,2)],'linestyle','--','color',[.6 .6 .6 .7],'linewidth',.2);
    errorbar([1,2],nanmean(SM),SSEM,'color',col,'linewidth',5);
    errorbar([1,2],nanmean(USM),USSEM,'color',col,'linestyle','--','linewidth',5);

ylabel('Log_{10} Mutual Info. [bits]');
set(gca,'xticklabel',{'Day 0','Day 1'});
set(gca,'xtick',[1:2],'tickdir','out','linewidth',2);
xlim([0.5,2.5]);
legend({'Stable','Unstable'},'Location','southwest');
if saveBool, print(linename,'-dpdf'); end

stability = [ones(nSessions*2,1); zeros(nSessions*2,1)];
session = repmat([ones(nSessions,1); 2*ones(nSessions,1)],[2,1]);
grps = {stability,session};
X = [SM(:);USM(:)];
[p,tbl,stats,terms] = anovan(X,grps,'model','full','varnames',{'Stability','Session'},...
    'display','off');
comps = multcompare(stats,'dimension',[1,2],'display','off','ctype','hsd');
disp(['Main effect of stability F(',num2str(tbl{2,3}),',',num2str(tbl{5,3}),...
        ') = ',num2str(tbl{2,6}),', P = ',num2str(tbl{2,7})]);
    disp(['Main effect of session F(',num2str(tbl{3,3}),',',num2str(tbl{5,3}),...
        ') = ',num2str(tbl{3,6}),', P = ',num2str(tbl{3,7})]);
    disp(['Interaction F(',num2str(tbl{4,3}),',',num2str(tbl{5,3}),...
        ') = ',num2str(tbl{4,6}),', P = ',num2str(tbl{4,7})]);
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