figure; hold on;
for i=215:2:219
    cd(MD(i).Location);
    load('TimeCells.mat', 'TimeCells')
    load('TemporalInfo.mat', 'sig')
    TimeCells = intersect(find(sig),TimeCells);
    
    corrStats = CorrTrdmllTrace(MD(i),MD(i+2),TimeCells);
    mapMD = getMapMD(MD(i));
    
    stable = msMatchCells(mapMD,MD(i:2:i+2),find(corrStats(:,2)<.01),false);
    unstable = msMatchCells(mapMD,MD(i:2:i+2),find(corrStats(:,2)>.01),false);
    
    [t1,m1] = getTimePeak(MD(i));
    [t2,m2] = getTimePeak(MD(i+2));

    h = scatter(m1(stable(:,1)),m2(stable(:,2)),10,'g','filled'); 
    alpha(h,.5);
    h = scatter(m1(unstable(:,1)),m2(unstable(:,2)),10,'k','filled');
    alpha(h,.5);
end

for i=253:255
    cd(MD(i).Location);
    load('TimeCells.mat', 'TimeCells')
    load('TemporalInfo.mat', 'sig')
    TimeCells = intersect(find(sig),TimeCells);
    
    corrStats = CorrTrdmllTrace(MD(i),MD(i+1),TimeCells);
    mapMD = getMapMD(MD(i));
    
    stable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)<.01),false);
    unstable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)>.01),false);
    
    [t1,m1] = getTimePeak(MD(i));
    [t2,m2] = getTimePeak(MD(i+1));

    h = scatter(m1(stable(:,1)),m2(stable(:,2)),10,'g','filled'); 
    alpha(h,.5);
    h = scatter(m1(unstable(:,1)),m2(unstable(:,2)),10,'k','filled');
    alpha(h,.5);
end

for i=274:277
    cd(MD(i).Location);
    load('TimeCells.mat', 'TimeCells')
    load('TemporalInfo.mat', 'sig')
    TimeCells = intersect(find(sig),TimeCells);
    
    corrStats = CorrTrdmllTrace(MD(i),MD(i+1),TimeCells);
    mapMD = getMapMD(MD(i));
    
    stable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)<.01),false);
    unstable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)>.01),false);
    
    [t1,m1] = getTimePeak(MD(i));
    [t2,m2] = getTimePeak(MD(i+1));

    h = scatter(m1(stable(:,1)),m2(stable(:,2)),10,'g','filled'); 
    alpha(h,.5);
    h = scatter(m1(unstable(:,1)),m2(unstable(:,2)),10,'k','filled');
    alpha(h,.5);
end

for i=287:290
    cd(MD(i).Location);
    load('TimeCells.mat', 'TimeCells')
    load('TemporalInfo.mat', 'sig')
    TimeCells = intersect(find(sig),TimeCells);
    
    corrStats = CorrTrdmllTrace(MD(i),MD(i+1),TimeCells);
    mapMD = getMapMD(MD(i));
    
    stable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)<.01),false);
    unstable = msMatchCells(mapMD,MD(i:i+1),find(corrStats(:,2)>.01),false);
    
    [t1,m1] = getTimePeak(MD(i));
    [t2,m2] = getTimePeak(MD(i+1));

    h = scatter(m1(stable(:,1)),m2(stable(:,2)),10,'g','filled'); 
    alpha(h,.5);
    h = scatter(m1(unstable(:,1)),m2(unstable(:,2)),10,'k','filled');
    alpha(h,.5);
end