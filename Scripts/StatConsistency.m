clear
loadMD;
fulldataset = MD(292:309); 
animals = unique({fulldataset.Animal});
nAnimals = length(animals);
colors = parula(nAnimals);

statType = 'si';
PCcrit = .01;

c = [];
stats = [];
for a=1:nAnimals
    N = 0;
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
    for s=1:length(ssns)-1
        cd(fulldataset(ssns(s)).Location);
        
        load('TimeCells.mat','TimeCells');
        load('TemporalInfo.mat','sig'); 
        
        %neurons = intersect(TimeCells,find(sig)); 
        neurons = getPlaceCells(fulldataset(ssns(s)),PCcrit);
        temp = msStats(fulldataset(ssns(s):ssns(s+1)),statType,neurons);
        
        stats = [stats; temp];
        
        N = N + size(temp,1); 
    end
    
    c = [c; repmat(colors(a,:),N,1)];
end

dpvalue = signrank(stats(:,1),stats(:,2));
[R,corrpvalue] = corr(stats(:,1),stats(:,2));

s = scatter(stats(:,1),stats(:,2),20,c,'filled');
alpha(s,.5);
set(gca,'tickdir','out');
axis tight;
yLims = get(gca,'ylim'); xLims = get(gca,'xlim');
%line(xLims,yLims,'color','k','linewidth',2);
l = lsline;
l.LineStyle = '--';
l.LineWidth = 2;
xlabel('Day 1'); ylabel('Day 2');
title('Mutual Spatial Information [bits]');
text(xLims(1),yLims(2),['p = ',num2str(corrpvalue)]);