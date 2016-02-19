function [MAP,MAPcols] = FilterMAPDates(batch_session_map,dates,sessionNums)
%[MAP,MAPcols] = FilterMAPDates(batch_session_map,dates,sessionNums)
%
%   This function will automatically scan batch_session_map.map for you
%   based on the other inputs. Its outputs are the column indices of MAP
%   that correspond to the dates and session numbers you indicated. 
%
%   INPUTS
%       batch_session_map: Output from neuron_reg_batch. 
%
%       dates: Cell array of dates in the form 'mm_dd_yyyy'. 
%
%       sessionNums: Vector, same size as dates indicating session number. 
%
%   OUTPUTS
%       MAP: batch_session_map.map, except with the first column removed.
%       
%       MAPcols: Column indices for MAP corresponding to the sessions
%       indicated in the inputs. 
%

%% Locate column indices. 
    %Get rid of the first column and NaNs. 
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