function MATRIX = PCEnsembleCorrMatrix(mds)
%
%
%

%%
    animals = unique({mds.Animal});     %Cell array of animal names. 
    nAnimals = length(animals); 
    nSessions = zeros(nAnimals,1);
    
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{mds.Animal})); 
        nSessions(a) = length(ssns);
    end
    nDays = max(nSessions);
    
    MATRIX = nan(nDays,nDays,nAnimals);
    for a=1:nAnimals
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal})); 
        nSessions = length(ssns);
        
        cellsOfInterest = cell(nSessions,1);
        for s=1:nSessions
            cd(mds(ssns(s)).Location);
            
            cellsOfInterest{s} = getPlaceCells(mds(ssns(s)),.01);
        end
        
        map = msMatchMultiSessionCells(mds(ssns),cellsOfInterest);
        
        for s1=1:nSessions
            neurons = map(:,s1);
            
            for s2=s1:nSessions
                neurons = EliminateUncertainMatches([mds(ssns(s1)) mds(ssns(s2))],neurons);
                
                temp = CorrPlaceFields(mds(ssns(s1)),mds(ssns(s2)),neurons);
                
                MATRIX(s1,s2,a) = nanmean(temp(:,1)); 
                MATRIX(s2,s1,a) = nanmean(temp(:,1)); 
            end
        end
                
    end
    
end