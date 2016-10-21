function matchMat = msMatchCells(mapMD,md,neurons)
%matchMat = msMatchCells(mapMD,md,neurons)
%
%   

%%
    load(fullfile(mapMD.Location,'batch_session_map.mat'));
    
    dates = {md.Date};
    sessions = [md.Session];
    nSessions = length(md);
    
    regDates = {batch_session_map.session.Date};
    regSessions = [batch_session_map.session.Session];
    
    MAP = batch_session_map.map(:,2:end);
    
    MAPcols = zeros(nSessions,1);
    for i=1:nSessions
        try
            MAPcols(i) = find(ismember(regDates,dates{i}) & ismember(regSessions,sessions(i)));
        catch
            error(['Error in above. Possible reason: MD input ', num2str(i),...
                ' has not been registered yet. Run neuron_reg_batch...']);
        end
    end
    
    neuronInd = ismember(MAP(:,MAPcols(1)),neurons);
    matchMat = MAP(neuronInd,MAPcols);
   
end
