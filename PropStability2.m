function [pct,stability,map,everStable] = PropStability2(mds,cellType)
%[pct,stability,map] = PropStability2(mds,cellType)
%
%

%% 
    nSessions = length(mds);
    
    %Load map.
%     mapMD = getMapMD(mds);
%     cd(mapMD.Location);
%     load('batch_session_map.mat');
%     map = batch_session_map.map(:,2:end); 
%     
%     %Reorder the map columns to reflect the order in mds. 
%     [~,~,mapCols] = msMatchCells(mapMD,mds,[],false);
%     map = map(:,mapCols);
    %Get all time/place cells. Only look at those that exist over two
    %consecutive days. 
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
        %Get critical p-value for stability. 
        temp = EliminateUncertainMatches([mds(s),mds(s+1)],temp);
        nNeuronsForStabilityCrit = length(temp);
        stableCrit = 0.01/nNeuronsForStabilityCrit;
        
        %Is a time/place cell on that day. 
        isCoding(:,s) = ismember(map(:,s),temp);
        
        %Do tuning curve correlations.
        switch cellType
            case 'time'
                corrs = CorrTrdmllTrace(mds(s),mds(s+1),neurons);
            case 'place'
                corrs = CorrPlaceFields(mds(s),mds(s+1),neurons);
        end
  
        %Find (un)stable cells.
        stable = find(corrs(:,2) < stableCrit);
        unstable = find(corrs(:,2) > stableCrit); 
        
        %Mark stable vs unstable. 
        i = ismember(map(:,s),stable);
        isStable(i,s) = 1;
        i = ismember(map(:,s),unstable);
        isStable(i,s) = 0;
                 
    end
    
    %Find out which cells on this day were time/place cells. 
    switch cellType
        case 'time'
            temp = getTimeCells(mds(end)); 
        case 'place'
            temp = getPlaceCells(mds(end),0.01);
    end
    isCoding(:,end) = ismember(map(:,end),temp); 
    
    %If cell was a time/place cell, but not stable, 0 it. 
    stableFlag = isCoding;
    stableFlag(isStable==0 & isCoding==1) = 0;
    
    %Get stable cells. 
    good = stableFlag(:,1:end-1) > 0;
    everStable = any(good,2);           %Cells that were ever stable. 
    %foundCell = map(:,1:end-1) > 0; 
    
    %If a cell was stable for all days it was detected, mark it as 1. 
    %pctDaysStable = sum(good,2) ./ sum(foundCell,2);
    nDaysDetected = sum(~isnan(isStable),2);
    pctDaysStable = sum(good,2) ./ nDaysDetected;
    stability.Stable = find(pctDaysStable == 1);

    %Get change in stability status. We define a cell as coming online on
    %the last day if it's a time/place cell on that day, but not on any
    %other days. 
    deltaStability = [zeros(size(stableFlag,1),1) diff(stableFlag,[],2)];
    deltaStability(:,end) = stableFlag(:,end) & all(stableFlag(:,1:end-1) == 0,2);
    stability.Incoming = find(any(ismember(deltaStability,1),2)); 
    
    %A cell is outgoing if its coding status goes from 1 to 0. 
    outgoing = find(any(ismember(deltaStability,-1),2)); 
    
    %Also handle ambiguous cases. That is, cells that don't change
    %stability status, but weren't classified as time/place cells. These
    %get labeled as outgoing because they were only time/place cells for
    %one day. 
    ambiguous = find(all(ismember(deltaStability,0),2));
    ambiguous = setdiff(ambiguous,stability.Stable);
    stability.Outgoing = union(outgoing,ambiguous);
    
    %We classified some incoming/outgoing cells as stable incorrectly.
    %Remove them from the stable list.
    ambiguous = intersect(union(stability.Outgoing,stability.Incoming),stability.Stable);
    stability.Stable(ismember(stability.Stable,ambiguous)) = [];
    
    %Get percentages. 
    pct.Stable = length(stability.Stable)/size(map,1);
    pct.Outgoing = length(stability.Outgoing)/size(map,1);
    pct.Incoming = length(stability.Incoming)/size(map,1);
end