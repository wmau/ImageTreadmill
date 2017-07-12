clear;
loadMD;

fulldataset = MD(292:309);
nSessions = length(fulldataset); 
cellType = 'placecells';
rasterType = 'place';
B = 500;

[skew,sig] = deal(cell(nSessions,1)); 
for s=1:nSessions
    disp(['Analyzing ',fulldataset(s).Animal,' ',fulldataset(s).Date]);
    
    skew{s} = getAllSkewnesses(fulldataset(s),'cellType',cellType,...
        'rasterType',rasterType); 
    sig{s} = nan(size(skew{s})); 
    neurons = AcquireTimePlaceCells(fulldataset(s),cellType)'; 
    
    shuffled = cell(B,1); 
    for i=1:B
        shuffled{i} = getAllSkewnesses(fulldataset(s),'cellType',cellType,...
            'rasterType',rasterType,'shuffle',true);
    end
    
    for n=neurons
        distribution = cellfun(@(x) x(n),shuffled); 
        
        p = sum(skew{s}(n) > distribution)/B; 
        
        if p < 0.025
            sig{s}(n) = 1;
        elseif p > 0.975
            sig{s}(n) = 2;
        else 
            sig{s}(n) = 0;
        end
    end
    

end

pEarly = sum(cellfun(@(x) nansum(x==1),sig))/sum(cell2mat(cellfun(@(x) sum(~isnan(x)),sig,'unif',0)));
pLate = sum(cellfun(@(x) nansum(x==2),sig))/sum(cell2mat(cellfun(@(x) sum(~isnan(x)),sig,'unif',0)));

disp(['Proportion of time cells that fire early: ',num2str(sum(cellfun(@(x) nansum(x==1),sig))),' ',num2str(pEarly*100),'%']);
disp(['Proportion of time cells that fire late: ',num2str(sum(cellfun(@(x) nansum(x==2),sig))),' ',num2str(pLate*100),'%']);