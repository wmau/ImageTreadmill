clear
loadMD;

fulldataset = MD(292:309); 
nSessions = length(fulldataset); 

[SI,TI] = deal([]);
for s=1:nSessions
    cd(fulldataset(s).Location); 
    load('Placefields_onTM.mat','RunOccMap','TMap_unsmoothed','exclude_frames');
    load('Pos_align.mat','PSAbool');
    PSAbool(:,exclude_frames) = false;
    
    TCs = AcquireTimePlaceCells(fulldataset(s),'dual'); 
    
    [I] = spatInfo(TMap_unsmoothed,RunOccMap,PSAbool,0); 
    
    SI = [SI; zscore(I(TCs))];
    
    load('TemporalInfo.mat','MI'); 
    TI = [TI; zscore(MI(TCs))];
end

scatter(TI,SI,'.');
[r,p] = corr(TI,SI);