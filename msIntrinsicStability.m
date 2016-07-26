function [pStableConnections,pStableCells] = msIntrinsicStability(mapMD,md,A)
%[pStableConnections,pStableCells] = msIntrinsicStability(mapMD,md,A)
%
%   Finds the proportion of targets and connections that are stable from
%   day to day. As of this writing, I use the KS-test MakeGraph method. 
%

%% Set up.
    nSessions = length(md);         %Number of sessions.
    [trigger,target] = find(A);     %Indices of triggers and targets for the initial session.
    nEdges = length(trigger);       %Number of connections.
    
    %Get the cell identities for each session.
    triggers = msMatchCells(mapMD,md,trigger); 
    targets = msMatchCells(mapMD,md,target);
    
    %Targets. 
    uTarg = unique(target);
    
    %Get the number of cells that mapped onto other sessions. 
    nGoodCells = zeros(1,nSessions); 
    for s=1:nSessions
        for e=uTarg'
            nGoodCells(s) = nGoodCells(s) + (targets(targets(:,1)==e,s) > 0);
        end
    end
    
    pStableConnections = zeros(1,nSessions); pStableConnections(1) = 1; 
    pStableCells = zeros(1,nSessions); pStableCells(1) = 1;
    goodCells = cell(1,nSessions);
    for s=2:nSessions
        load(fullfile(md(s).Location,'graphData_KS.mat'),'A');
        n = 0; nGoodConnections = 0;
        for e=1:nEdges
            thisTrig = trigger(e);
            thisTarg = target(e); 
            
            one = triggers(triggers(:,1)==thisTrig,s);
            two = targets(targets(:,1)==thisTarg,s);
            
            if one > 0 && two > 0
                n = n+A(one,two);
                nGoodConnections = nGoodConnections + 1;
                
                if A(one,two)
                    goodCells{s} = [goodCells{s} thisTarg];
                end
            end
        end
        
        pStableConnections(s) = n/nGoodConnections;
        pStableCells(s) = length(unique(goodCells{s}))/nGoodCells(s);
    end

    
end