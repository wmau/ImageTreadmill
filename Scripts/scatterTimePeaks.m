loadMD;
fulldataset = [MD(215:2:221) MD(253:256) MD(274:278) MD(287:291)];
animals = unique({fulldataset.Animal});

[~,~,stable,unstable] = PartitionStats(fulldataset,'time','TI');
nAnimals = length(stable);

figure; hold on;
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions = length(ssns);
    mapMD = getMapMD(fulldataset(ssns(1)));
    
    for s=1:nSessions-1
        stableMatches = msMatchCells(mapMD,fulldataset(ssns(s:s+1)),stable{a}{s},true);

        unstableMatches = msMatchCells(mapMD,fulldataset(ssns(s:s+1)),unstable{a}{s},true);
        
        [~,t1] = getTimePeak(fulldataset(ssns(s)));
        [~,t2] = getTimePeak(fulldataset(ssns(s+1)));
        
        scatter(t1(stableMatches(:,1)),t2(stableMatches(:,2)),10,...
            [0 .5 .5],'filled');
        scatter(t1(unstableMatches(:,1)),t2(unstableMatches(:,2)),10,...
            [.7 .7 .7],'filled');
    end
    
end