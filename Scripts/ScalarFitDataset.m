clear;
clc;

loadMD;
fulldataset = MD(292:309);
nSessions = length(fulldataset); 
animals = unique({fulldataset.Animal});
nAnimals = length(animals);

for s=1:nSessions
    PC_DATA(s).peaks = getPlaceCellPeak(fulldataset(s));
    PC_DATA(s).Date = fulldataset(s).Date;
    PC_DATA(s).Animal = fulldataset(s).Animal;
    PC_DATA(s).PlaceCells = getPlaceCells(fulldataset(s),0.01);
end

map = cell(nAnimals,1);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    nSessions = length(ssns); 
    
    neurons = cell(nSessions,1);
    for s=1:nSessions
        neurons{s} = getPlaceCells(fulldataset(ssns(s)),0.01); 
    end
    
    map{a} = msMatchMultiSessionCells(fulldataset(ssns),neurons);
end