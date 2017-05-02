function PropStability(mds,cellType)
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
    end
    map = msMatchMultiSessionCells(mds,cellsOfInterest);

%%
    %Preallocate stability matrix. 
    isStable = nan(size(map));
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
    keyboard;
end