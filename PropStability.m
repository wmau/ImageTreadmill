function [pct,stability,map] = PropStability(mds,cellType)
%
%
%

%% 
    nSessions = length(mds);
    
    %Load map.
    mapMD = getMapMD(mds);
%     cd(mapMD.Location);
%     load('batch_session_map.mat');
%     map = batch_session_map.map(:,2:end); 
%     
%     %Reorder the map columns to reflect the order in mds. 
%     [~,~,mapCols] = msMatchCells(mapMD,mds,[],false);
%     map = map(:,mapCols);
    cellsOfInterest = cell(nSessions,1);
    for s=1:nSessions
        switch cellType
            case 'time'
                cellsOfInterest{s} = getTimeCells(mds(s));
            case 'place'
                cellsOfInterest{s} = getPlaceCells(mds(s),0.01);
        end
        if s < nSessions
            cellsOfInterest{s} = EliminateUncertainMatches([mds(s),mds(s+1)],cellsOfInterest{s});
        end
    end
    map = msMatchMultiSessionCells(mds,cellsOfInterest);

%%
    %Preallocate stability matrix. 
    [isCoding,isStable] = deal(nan(size(map)));
    for s=1:nSessions-1
        cd(mds(s).Location); 
        
        %Get neuron IDs from this session. Only look at cells that map to
        %the next day. 
        neurons = map(map(:,s+1)>0,s);
        neurons(neurons==0) = [];
                 
        %Get cell numbers for stability criterion. 
        switch cellType
            case 'time'
                temp = getTimeCells(mds(s));
            case 'place'
                temp = getPlaceCells(mds(s),0.01);
        end
        temp = EliminateUncertainMatches([mds(s),mds(s+1)],temp);
        nNeuronsForStabilityCrit = length(temp);
        stableCrit = 0.01/nNeuronsForStabilityCrit;
        
        isCoding(:,s) = ismember(map(:,s),temp);
        
        %Do tuning curve correlations.
        switch cellType
            case 'time'
                corrs = CorrTrdmllTrace(mds(s),mds(s+1),neurons);
            case 'place'
                corrs = CorrPlaceFields(mds(s),mds(s+1),neurons);
        end
  
        stable = find(corrs(:,2) < stableCrit);
        unstable = find(corrs(:,2) > stableCrit); 
        
        i = ismember(map(:,s),stable);
        isStable(i,s) = 1;
        i = ismember(map(:,s),unstable);
        isStable(i,s) = 0;
                 
    end
    
    %Fill in the last row for isCoding
    switch cellType
        case 'time'
            temp = getTimeCells(mds(end)); 
        case 'place'
            temp = getPlaceCells(mds(end),0.01);
    end
    isCoding(:,end) = ismember(map(:,end),temp); 
    
    good = isStable(:,1:end-1) > 0;
    %foundCell = map(:,1:end-1) > 0; 
    
    %pctDaysStable = sum(good,2) ./ sum(foundCell,2);
    pctDaysStable = sum(good,2) ./ (size(isStable,2)-1);
    stability.Stable = find(pctDaysStable == 1);
    
    deltaStability = [zeros(size(isStable,1),1) diff(isStable,[],2)];
    deltaStability(:,end) = isCoding(:,end) & all(isCoding(:,1:end-1) == 0,2);
    deltaStability(isnan(deltaStability)) = 0;
    stability.Incoming = find(any(ismember(deltaStability,1),2)); 
    outgoing = find(any(ismember(deltaStability,-1),2)); 
    
    ambiguous = find(all(ismember(deltaStability,0) | isnan(deltaStability),2));
    ambiguous = setdiff(ambiguous,stability.Stable);
    stability.Outgoing = union(outgoing,ambiguous);
    
    ambiguous = intersect(union(stability.Outgoing,stability.Incoming),stability.Stable);
    stability.Stable(ismember(stability.Stable,ambiguous)) = [];
    
    pct.Stable = length(stability.Stable)/size(map,1);
    pct.Outgoing = length(stability.Outgoing)/size(map,1);
    pct.Incoming = length(stability.Incoming)/size(map,1);
    
end