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
    cd(mds(1).Location); 
    load('FinalOutput.mat','NumNeurons');

    %Match cells. 
    matches = msMatchCells(mds,neurons,true);
    
    stats = nan(NumNeurons,nSessions);
    %Toss into matrix. 
    for s=1:nSessions
        stats(matches(:,1),s) = DATA.(stat){s}(matches(:,s));
    end
    %stats(matches(:,s),s) = DATA.
    
end