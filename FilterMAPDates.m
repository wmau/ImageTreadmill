function [MAP,MAPcols] = FilterMAPDates(batch_session_map,dates,sessionNums)
%
%
%

%%
    MAP = batch_session_map.map(:,2:end); MAP(isnan(MAP)) = 0; 
    
    %Concatenate dates and session numbers in the registered sessions. 
    MAPdates = {batch_session_map.session.Date}; 
    MAPsessionnums = [batch_session_map.session.Session];
    nSessions = length(sessionNums);
    
    %Find columns that correspond to the dates in base and comp. 
    MAPcols = zeros(nSessions,1); 
    for i=1:nSessions
        try
            MAPcols(i) = find(ismember(MAPdates,dates{i}) ...
                & ismember(MAPsessionnums,sessionNums(i)));
        catch
            error(['Error in above. Possible reason: MD input has not been registered yet. '...
                'Run neuron_reg_batch...']);
        end
    end
    
end