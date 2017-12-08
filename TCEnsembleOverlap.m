function [TCsonBothDays,nTCsInBothSessions,p] = ...
    TCEnsembleOverlap(base,comp,reference,rSamp)
%[TCsonBothDays,nTCsInBothSessions,p] = ...
%    TCEnsembleOverlap(base,comp,reference)
%
%   

%%  
    baseTCs = getTimeCells(base); 
    
    if rSamp
        load(fullfile(comp.Location,'FinalOutput.mat'),'NumNeurons');
        nTCs = length(getTimeCells(comp)); 
        compTCs = randsample(NumNeurons,nTCs);
    else 
        compTCs = getTimeCells(comp); 
    end
    
    matchedCells = msGetMap([base,comp]); 
    
    rowsDay1TCs = ismember(matchedCells(:,1),baseTCs);  %Rows of the map that were TCs on day 1.
    baseTCsonDay2 = matchedCells(rowsDay1TCs,2);        %Cell #s of day 1 TCs on day 2. 
    day1and2TCs = ismember(baseTCsonDay2,compTCs);      %Indices of day 1 TCs that were also TCs on day 2. 
    
    %Cells that were TCs on both days. Columns correspond to different
    %sessions.
    TCsonBothDays = matchedCells(ismember(matchedCells(:,2),...
        baseTCsonDay2(day1and2TCs)),:);

    %Number of these cells. 
    nTCsInBothSessions = length(TCsonBothDays); 
    
    %To calculate the percentage of overlap, the denominator is ambiguous.
    %Here are some possible choices. They are listed in order of largest to
    %smallest pool. 
    TCMap = msMatchMultiSessionCells([base,comp],{baseTCs compTCs});
    switch reference
        case 'allTCs'                   %ALL time cells, regardless of whether they were found on both or one day.
            n = size(TCMap,1); 
        case 'allMappedTCs'             %ALL time cells that also exist across both days. 
            n = sum(all(TCMap>0,2));
        case 'day1TCs'                  %Time cells on day 1, regardless of whether they were found on day 2. 
            day1TCMap = msMatchCells([base,comp],baseTCs,false);
            
            n = size(day1TCMap,1); 
        case 'day1MappedTCs'            %Time cells on day 1 that were also found on day 2. 
            day1TCMap = msMatchCells([base,comp],baseTCs,true); 
            
            n = size(day1TCMap,1);
    end
    p = nTCsInBothSessions./n;
end