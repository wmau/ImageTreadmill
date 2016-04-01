function pMatched = nMatchedNeurons(mapMD,MD1,MDs)
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
    pMatched = sum(MAP(found,MAPcols(2:end))>0)./length(found);
    
end