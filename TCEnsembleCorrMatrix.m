function TCEnsembleCorrMatrix(mds)
%
%
%

%%
    animals = unique({mds.Animal});     %Cell array of animal names. 
    nAnimals = length(animals); 
    
    for a=1:nAnimals
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal})); 
        nSessions = length(ssns);
        
        
    end
end