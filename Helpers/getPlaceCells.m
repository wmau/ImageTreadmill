function PCs = getPlaceCells(md,crit)
%PCs = getPlaceCells(md,crit)
%
%   Gets all the cells with significant mutual spatial information. 
%
%   INPUTS
%       md: Session entry.
%
%       crit: p-value for consideration as a place cell.
%
%   OUTPUT
%       PCs: Indices of place cells. 

%% Get place cells.
    cd(md.Location);
    load('Placefields.mat','TMap_gauss','pval');
    load('PlacefieldStats.mat','PFnHits','bestPF');
    load('SpatialInfo.mat','MI');
    
    idx = sub2ind(size(PFnHits),1:size(PFnHits,1),bestPF');    
    PCs = find(pval<crit & MI'>0 & PFnHits(idx) > 10);
    
end