function [nMatched,pMatched,found,MAP,MAPcols] = nMatchedNeurons(mapMD,MD1,MDs)
%
%
%

%% 
    s = [MD1,MDs];
    dates = {s.Date};
    sNums = [s.Session];
    
    load(fullfile(mapMD.Location,'batch_session_map.mat'));
    [MAP,MAPcols] = FilterMAPDates(batch_session_map,dates,sNums); 
   
    found = find(MAP(:,MAPcols(1)));
    nMatched = sum(MAP(found,MAPcols(2:end))>0);
    pMatched = nMatched./length(found);
    
end