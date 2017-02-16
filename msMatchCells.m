function [matchMat,mapRows,mapCols] = msMatchCells(mapMD,md,neurons,trim)
%matchMat = msMatchCells(mapMD,md,neurons)
%
%   Matches all specified neurons across sessions specified in md. 
%
%   INPUTS
%       mapMD: MD entry containing batch_session_map.mat.
%
%       md: Specified sessions.
%
%       neurons: Specified neurons.
%
%       trim: Logical, whether or not you want to excise unmapped neurons.
%
%   OUTPUT
%       matchMat: NxS matrix (N = # neurons, S = # sessions). First column
%       is the neurons vector. Subsequent columns are the corresponding
%       neurons from other sessions.
%

%% Match cells.
    load(fullfile(mapMD.Location,'batch_session_map.mat'));
    
    %Specified sessions.
    dates = {md.Date};
    sessions = [md.Session];
    nSessions = length(md);
    
    %Sessions registered.
    regDates = {batch_session_map.session.Date};
    regSessions = [batch_session_map.session.Session];
    
    %Trim first column.
    MAP = batch_session_map.map(:,2:end);
    
    %Find columns in the mapping matrix that corresponding to the specified
    %sessions.
    mapCols = zeros(nSessions,1);
    for i=1:nSessions
        try
            mapCols(i) = find(ismember(regDates,dates{i}) & ismember(regSessions,sessions(i)));
        catch
            error(['Error in above. Possible reason: MD input ', num2str(i),...
                ' has not been registered yet. Run neuron_reg_batch...']);
        end
    end
    
    %Get the corresponding rows.
    mapRows = find(ismember(MAP(:,mapCols(1)),neurons));
    matchMat = MAP(mapRows,mapCols);
   
    if trim
        matchMat(any(matchMat(:,2:end)==0,2),:) = [];
        matchMat(any(isnan(matchMat(:,2:end)),2),:) = [];
    end
end
