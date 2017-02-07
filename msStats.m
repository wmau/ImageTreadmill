function stats = msStats(mds,stat,neurons)
%stats = msStats(mds,stat,neurons)
%
%   Gets an arbitrary statistic and gathers it from each session in mds.
%   Only tested for spatial and temporal information, so these statistics
%   should be scalars. 
%
%   INPUTS
%       mds: Sessions you are looking at. 
%
%       stat: String, must match an acceptable field in
%       CompileMultiSessionData. 
%
%       neurons: Indices for the first session in mds. 
%
%   OUTPUT
%       stats: NxS matrix (N=number of mapped neurons, S=number of sessions
%       in mds) containing all your statistics. 
%   

%% Compile statistics. 
    stat = lower(stat);
    DATA = CompileMultiSessionData(mds,{stat});
    nSessions = length(mds);
    nNeurons = length(neurons);
    
    %Session containing mapping matrix. 
    mapMD = getMapMD(mds);
    
    %Match cells. 
    matches = msMatchCells(mapMD,mds,neurons,true);
    
    stats = nan(nNeurons,nSessions);
    %Toss into matrix. 
    stats(matches(:,1),1) = DATA.(stat){1}(matches(:,1));
    stats(matches(:,1),2) = DATA.(stat){2}(matches(:,2));
    %stats(matches(:,s),s) = DATA.
    
end