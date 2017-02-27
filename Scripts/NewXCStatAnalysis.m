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

[new1,new2,r1,r2,N1,N2,R1,R2] = deal([]);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions = length(ssns)-1;
    
    for s=1:nSessions
        
        switch cellType
            case 'Time'
                [stat1,random1,stat2,random2] = NewTimeCellStats(fulldataset(ssns(s)),...
                    fulldataset(ssns(s+1)),statType,'plotit',false);
            case 'Place'
                [stat1,random1,stat2,random2] = NewPlaceCellStats(fulldataset(ssns(s)),...
                    fulldataset(ssns(s+1)),statType,'plotit',false);
        end
        
        new1 = [new1; stat1];
        new2 = [new2; stat2];
        r1 = [r1; random1];
        r2 = [r2; random2];
        
        N1 = [N1; mean(new1)];
        N2 = [N2; mean(new2)];
        R1 = [R1; mean(r1)];
        R2 = [R2; mean(r2)];
    end
end

%Histogram of statistic. 
figure; hold on
histogram(new1,'normalization','probability','binwidth',binsize,'edgecolor','none','facecolor',c);
histogram(r1,'normalization','probability','binwidth',binsize,'edgecolor','none','facecolor',[.7 .7 .7]);
xlabel(yLabel); 
ylabel('Proportion');
legend({['Future ',cellType, ' Cells'],'Random Sample'});

figure('Position',[440 350 250 420]); hold on;
m = [mean(N1) mean(N2)];
Rm = [mean(R1) mean(R2)];
sem = [std(N1)./sqrt(length(N1)) std(N2)./sqrt(length(N2))];
bar(1,m(1),'facecolor',c,'edgecolor',sc,'linewidth',3,'facealpha',.5);
b = bar(2,m(2),'facecolor',c,'edgecolor',sc,'linewidth',3,'facealpha',.9);
w = b.BarWidth;
for i=1:2
    line([i-w/2 i+w/2],[Rm(i) Rm(i)],'color','b','linewidth',2);
    errorbar(i,m(i),sem(i),'linewidth',2,'color',sc);
end

for i=1:length(N1)
    plot([1,2],[N1(i) N2(i)],'color',[sc .5]);
end

set(gca,'xticklabel',{'Day -1','Day 0'});
set(gca,'xtick',[1:2],'tickdir','out','linewidth',2);
ylabel(yLabel);
xlim([0.5,2.5]);
pN1Random = ranksum(N1,R1);
pN2Random = ranksum(N2,R2);
pN1N2 = signrank(N1,N2); 
title({['Pre-Crit vs R, P = ',num2str(pN1Random)],...
       ['Pre-Crit vs Crit, P = ',num2str(pN1N2)],...
       ['Crit vs R, P = ',num2str(pN2Random)]});
if saveBool
    print(dayComp,'-dpdf');
end
