function stats = corrStats(mds,stat1,stat2)
%stats = corrStats(mds,stat1,stat2)
%
%   Gets two different types of stats for each specified session and just
%   makes a long two-column matrix where column 1 is stat1 and column 2 is
%   stat2. Currently only looks at stats(stats>0) for some reason I don't
%   remember anymore. Should probably be adapted to normalize within
%   session.
%

%% 
    DATA = CompileMultiSessionData(mds,{stat1,stat2});
    nSessions = length(mds);
    
    stats = [];
    for s=1:nSessions
        stats = [   stats; 
                    DATA.(stat1){s}, DATA.(stat2){s}];
    end
    
    good = stats(:,1)>0 & stats(:,2)>0;
    scatter(stats(good,1),stats(good,2),'.');
    keyboard; 
end