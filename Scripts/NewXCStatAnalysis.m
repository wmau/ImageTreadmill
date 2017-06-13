clear
loadMD
%MD(292:295) = G45.
%MD(296:299) = G48.
%MD(300:304) = Bellatrix.
%MD(305:309) = Polaris.
fulldataset = MD(292:309);   
animals = unique({fulldataset.Animal});
nAnimals = length(animals);

statType = 'si';
cellType = 'Time';

saveBool = true;
folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures';
dayComp = fullfile(folder,['New ',cellType,' Cell ',statType,' Day 1 vs Day 2']);

if saveBool
    c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');
    
    if ~strcmp(c,'y')
        saveBool = false;
    end
end

binsize = .001;
switch statType
    case 'ti', yLabel = 'Temporal Info.';
    case 'si', yLabel = 'Spatial Info.';
    case 'fr', yLabel = 'Ca Event Freq.';
end 

switch cellType
    case 'Place', c = [.58 .44 .86];
    case 'Time', c = [0 .5 .5];
end

switch statType
    case 'ti', sc = [0 .5 .5];
    case 'si', sc = [.58 .44 .86]; 
    case 'fr', sc = [.7 .7 .7];
end

rc = [.5 .5 .5];

[newStat1,newStat2,popStat1,popStat2,r1,r2,...
    newMean1,newMean2,popMean1,popMean2,R1,R2] = deal([]);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions = length(ssns)-1;
    
    for s=1:nSessions
        
        switch cellType
            case 'Time'
                [newStatS1,random1,stat1,newStatS2,random2,stat2] = NewTimeCellStats(fulldataset(ssns(s)),...
                    fulldataset(ssns(s+1)),statType,'plotit',false);
            case 'Place'
                [newStatS1,random1,stat1,newStatS2,random2,stat2] = NewPlaceCellStats(fulldataset(ssns(s)),...
                    fulldataset(ssns(s+1)),statType,'plotit',false);
        end
        
        newStat1 = [newStat1; newStatS1];
        newStat2 = [newStat2; newStatS2];
        
        popStat1 = [popStat1; stat1];
        popStat2 = [popStat2; stat2];
        
        r1 = [r1; random1];
        r2 = [r2; random2];
        
        
        newMean1 = [newMean1; mean(newStatS1)];
        newMean2 = [newMean2; mean(newStatS2)];
        
        popMean1 = [popMean1; mean(stat1)];
        popMean2 = [popMean2; mean(stat2)];
        
        R1 = [R1; mean(random1)];
        R2 = [R2; mean(random2)];
    end
end

newStat1 = log10(newStat1);
newStat2 = log10(newStat2);
popStat1 = log10(popStat1); 
popStat2 = log10(popStat2);
r1 = log10(r1);
r2 = log10(r2);
newMean1 = log10(newMean1); 
newMean2 = log10(newMean2);
R1 = log10(R1);
R2 = log10(R2);

%Histogram of statistic. 
figure; hold on
histogram(newStat1,'normalization','probability','binwidth',binsize,'edgecolor','none','facecolor',c);
histogram(r1,'normalization','probability','binwidth',binsize,'edgecolor','none','facecolor',[.7 .7 .7]);
xlabel(yLabel); 
ylabel('Proportion');
legend({['Future ',cellType, ' Cells'],'Random Sample'});

figure('Position',[440 350 250 420]); hold on;
m = [mean(newMean1) mean(newMean2)];
Rm = [mean(R1) mean(R2)];
sem = [std(newMean1)./sqrt(length(newMean1)) std(newMean2)./sqrt(length(newMean2))];
bar(1,m(1),'facecolor',c,'edgecolor',sc,'linewidth',3,'facealpha',.5);
b = bar(2,m(2),'facecolor',c,'edgecolor',sc,'linewidth',3,'facealpha',.9);
w = b.BarWidth;
for i=1:2
    line([i-w/2 i+w/2],[Rm(i) Rm(i)],'color','b','linewidth',2);
    errorbar(i,m(i),sem(i),'linewidth',2,'color',sc);
end

for i=1:length(newMean1)
    plot([1,2],[newMean1(i) newMean2(i)],'color',[sc .5]);
end

set(gca,'xticklabel',{'Day -1','Day 0'});
set(gca,'xtick',[1:2],'tickdir','out','linewidth',4);
ylabel(yLabel);
xlim([0.5,2.5]);
pN1Random = ranksum(newMean1,R1);
pN2Random = ranksum(newMean2,R2);
pN1N2 = signrank(newMean1,newMean2); 
title({['Pre-Crit vs R, P = ',num2str(pN1Random)],...
       ['Pre-Crit vs Crit, P = ',num2str(pN1N2)],...
       ['Crit vs R, P = ',num2str(pN2Random)]});


%%
nSessions = length(newMean1);
newMeans = [nanmean(newMean1) nanmean(newMean2)];
newSEMs = [nanstd(newMean1)/sqrt(nSessions) nanstd(newMean2)/sqrt(nSessions)];

randMeans = [nanmean(R1) nanmean(R2)];
randSEMs = [nanstd(R1)/sqrt(nSessions) nanstd(R2)/sqrt(nSessions)];

figure('Position',[830 240 250 440]); hold on;
for i=1:length(newMean1)
    plot([1,2],[newMean1(i) newMean2(i)],'color',[sc .5]);
    plot([1,2],[R1(i) R2(i)],'color',[0 0 1 .5]);
end
errorbar([1,2],newMeans,newSEMs,'color',sc,'linewidth',5);
errorbar([1,2],randMeans,randSEMs,'color','b','linewidth',5);
xlim([0.5,2.5]);
set(gca,'xticklabel',{'Day -1','Day 0'},'xtick',[1:2],'tickdir','out','linewidth',4,'fontsize',15);
ylabel(yLabel);
if saveBool
    print(dayComp,'-dpdf');
end

%%
X = [newMean1; newMean2; R1; R2];
time = repmat([ones(nSessions,1); 2*ones(nSessions,1)],[2,1]);
sample = [ones(nSessions*2,1); 2*ones(nSessions*2,1)];
grps = {time,sample};
[~,tbl,stats] = anovan(X,grps,'model','full','varnames',{'Day','Sample'},'display','off');
comps = multcompare(stats,'dimension',[1,2],'display','off');
 disp(['Main effect of day F(',num2str(tbl{2,3}),',',num2str(tbl{5,3}),...
        ') = ',num2str(tbl{2,6}),', P = ',num2str(tbl{2,7})]);
    disp(['Main effect of sample F(',num2str(tbl{3,3}),',',num2str(tbl{5,3}),...
        ') = ',num2str(tbl{3,6}),', P = ',num2str(tbl{3,7})]);
    disp(['Interaction F(',num2str(tbl{4,3}),',',num2str(tbl{5,3}),...
        ') = ',num2str(tbl{4,6}),', P = ',num2str(tbl{4,7})]);
    disp(['Difference between Day -1 and Day 0 in new coding cells P = ',num2str(comps(1,6))]);
    disp(['Difference between Day -1 and Day 0 in random cells P = ',num2str(comps(6,6))]);
    disp(['Difference between new coding cells and random on Day -1 P = ',num2str(comps(2,6))]);
    %disp(['Difference between new coding cells and random on Day 0 P = ',num2str(comps(5,6))]);