function [ALLONSETS,plotme] = triggerRanks(md)
%[ALLONSETS,plotme] = triggerRanks(md)
%
%   Plots 

%%
    DATA = CompileMultiSessionData(md,{'t','timecells','ft','ttl'}); 
    T = unique(cell2mat(DATA.t));
    assert(length(T)==1,...
        'Error: Variable treadmill run durations present in this dataset.');
    
    nSessions = length(md); 
    animals = unique({md.Animal});
    nAnimals = length(unique({md.Animal})); 
    plotme = zeros(T,nAnimals); 
    %allpeaks = [];
    edges = linspace(0,10,11);
    ALLONSETS = [];
    
    for thisSession = 1:nSessions
        inds = DATA.ttl{thisSession}.inds(logical(DATA.ttl{thisSession}.complete),:);
        inds(:,2) = inds(:,1) + 20*DATA.t{thisSession}-1; 
        
        cd(md(thisSession).Location); 
        load('graphData_p','A'); 
        
        [trigger,target] = find(A); 
        nConnections = length(trigger);
        TMAlignedOnsets = [];
        for i=1:nConnections
            immRaster = buildRaster(inds,DATA.ft{thisSession},trigger(i));
            targRaster = buildRaster(inds,DATA.ft{thisSession},target(i)); 
            
            onsets = TMLatencies(immRaster,targRaster); 
            
            TMAlignedOnsets = [TMAlignedOnsets, onsets];
            ALLONSETS = [ALLONSETS, onsets]; 
        end
        %target = unique(target); 
        %target = intersect(target,TCdata.timecells{thisSession});
        
%         tPeaks = getTimePeak(md(thisSession)); 
%         allpeaks = [allpeaks; tPeaks(target)];
%         
%         count = histcounts(tPeaks(target),linspace(0,10,10));

        count = histcounts(TMAlignedOnsets,edges);
        
        animalind = find(strcmp(animals,md(thisSession).Animal));      
        try
            plotme(:,animalind) = plotme(:,animalind) + count;
        catch
            plotme(:,animalind) = plotme(:,animalind) + count'; 
        end
    end
    
    figure;
    bar(plotme,'stacked');
end