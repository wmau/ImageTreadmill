clear
loadMD
%MD(292:295) = G45.
%MD(296:299) = G48.
%MD(300:304) = Bellatrix.
%MD(305:309) = Polaris.
fulldataset = MD(292:309);   
animals = unique({fulldataset.Animal});
nAnimals = length(animals);

statType = 'fr';
cellType = 'Place';

binsize = .1;
switch statType
    case 'ti', yLabel = 'Mean Norm. Temporal Info. [z-bits]';
    case 'si', yLabel = 'Mean Norm. Spatial Info. [z-bits]';
    case 'fr', yLabel = 'Mean Norm. Ca Event Freq.';
end

switch cellType
    case 'Place', c = [.58 .44 .86];
    case 'Time', c = [0 .5 .5];
end

[new,r,N,R] = deal([]);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions = length(ssns)-1;
    
    for s=1:nSessions
        
        switch cellType
            case 'Time'
                [stat,randomSample] = NewTimeCellStats(fulldataset(ssns(s)),...
                    fulldataset(ssns(s+1)),statType,'plotit',false);
            case 'Place'
                [stat,randomSample] = NewPlaceCellStats(fulldataset(ssns(s)),...
                    fulldataset(ssns(s+1)),statType,'plotit',false);
        end
        
        new = [new; stat];
        r = [r; randomSample];
        
        N = [N; mean(new)];
        R = [R; mean(r)];
    end
end

%Histogram of statistic. 
figure; hold on
histogram(new,'normalization','probability','binwidth',binsize,'edgecolor','none','facecolor',c);
histogram(r,'normalization','probability','binwidth',binsize,'edgecolor','none','facecolor',[.7 .7 .7]);
xlabel(yLabel(6:end)); 
ylabel('Proportion');
legend({['Future ',cellType, ' Cells'],'Random Sample'});

figure('Position',[400 240 250 440]); hold on;
m = [mean(N) mean(R)];
sem = [std(N)./sqrt(length(N)) std(R)./sqrt(length(R))];
bar(1,m(1),'facecolor',c,'edgecolor','none','linewidth',3);
bar(2,m(2),'facecolor',[.7 .7 .7],'edgecolor','none','linewidth',3);
for i=1:length(N)
    plot([1 2],[N(i) R(i)],'color',[0 0 0 .3],'linewidth',.5);
end
errorbar(1,m(1),sem(1),'linewidth',2,'color',c);
errorbar(2,m(2),sem(2),'linewidth',2,'color',[.7 .7 .7]);
set(gca,'xticklabel',{['Future ',cellType,' Cell'],'Random'});
set(gca,'xtick',[1:2],'tickdir','out','linewidth',2);
ylabel(yLabel);
xlim([0.5,2.5]);
p = signrank(N,R); 
title(['P = ',num2str(p)]);