figure; hold on;
for i=292:294
    cd(MD(i).Location);
    load('TimeCells.mat', 'TimeCells')
    load('TemporalInfo.mat', 'sig')
    TimeCells = intersect(find(sig),TimeCells);
    
    corrStats = CorrTrdmllTrace(MD(i),MD(i+1),TimeCells);
    mapMD = getMapMD(MD(i));
    
    crit = .01/length(TimeCells);
    
    stable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)<crit),false);
    unstable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)>crit),false);
    
    [t1,m1] = getTimePeak(MD(i));
    [t2,m2] = getTimePeak(MD(i+1));

    h = scatter(m1(stable(:,1)),m2(stable(:,2)),10,'g','filled'); 
    alpha(h,.5);
    h = scatter(m1(unstable(:,1)),m2(unstable(:,2)),10,'k','filled');
    alpha(h,.5);
end

for i=296:298
    cd(MD(i).Location);
    load('TimeCells.mat', 'TimeCells')
    load('TemporalInfo.mat', 'sig')
    TimeCells = intersect(find(sig),TimeCells);
    
    corrStats = CorrTrdmllTrace(MD(i),MD(i+1),TimeCells);
    mapMD = getMapMD(MD(i));
    
    crit = .01/length(TimeCells);
    
    stable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)<crit),false);
    unstable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)>crit),false);
    
    [t1,m1] = getTimePeak(MD(i));
    [t2,m2] = getTimePeak(MD(i+1));

    h = scatter(m1(stable(:,1)),m2(stable(:,2)),10,'g','filled'); 
    alpha(h,.5);
    h = scatter(m1(unstable(:,1)),m2(unstable(:,2)),10,'k','filled');
    alpha(h,.5);
end

for i=300:303
    cd(MD(i).Location);
    load('TimeCells.mat', 'TimeCells')
    load('TemporalInfo.mat', 'sig')
    TimeCells = intersect(find(sig),TimeCells);
    
    corrStats = CorrTrdmllTrace(MD(i),MD(i+1),TimeCells);
    mapMD = getMapMD(MD(i));
    
    crit = .01/length(TimeCells);
    
    stable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)<crit),false);
    unstable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)>crit),false);
    
    [t1,m1] = getTimePeak(MD(i));
    [t2,m2] = getTimePeak(MD(i+1));

    h = scatter(m1(stable(:,1)),m2(stable(:,2)),10,'g','filled'); 
    alpha(h,.5);
    h = scatter(m1(unstable(:,1)),m2(unstable(:,2)),10,'k','filled');
    alpha(h,.5);
end

for i=305:308
    cd(MD(i).Location);
    load('TimeCells.mat', 'TimeCells')
    load('TemporalInfo.mat', 'sig')
    TimeCells = intersect(find(sig),TimeCells);
    
    corrStats = CorrTrdmllTrace(MD(i),MD(i+1),TimeCells);
    mapMD = getMapMD(MD(i));
    
    crit = .01/length(TimeCells);
    
    stable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)<crit),false);
    unstable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)>crit),false);
    
    [t1,m1] = getTimePeak(MD(i));
    [t2,m2] = getTimePeak(MD(i+1));

    h = scatter(m1(stable(:,1)),m2(stable(:,2)),10,'g','filled'); 
    alpha(h,.5);
    h = scatter(m1(unstable(:,1)),m2(unstable(:,2)),10,'k','filled');
    alpha(h,.5);
end
line([0 10], [0 10],'color','k','linestyle',':');
xlabel('Time Peak, Day 1 [s]'); ylabel('Time Peak, Day 2 [s]');
set(gca,'tickdir','out');