function msIntrinsicStability(mapMD,md,A)
%
%
%

%%
    nSessions = length(md);
    [trigger,target] = find(A);
    nEdges = length(trigger); 
    
    triggers = msMatchCells(mapMD,md,trigger); 
    targets = msMatchCells(mapMD,md,target);
    
    uTarg = unique(target);
    
    nGoodCells = zeros(1,nSessions); 
    for s=1:nSessions
        for e=uTarg'
            try 
            nGoodCells(s) = nGoodCells(s) + (targets(targets(:,1)==e,s) > 0);
            catch
                keyboard;
            end
        end
    end
    
    pStableConnections = zeros(1,nSessions); pStableConnections(1) = 1; 
    pStableCells = zeros(1,nSessions); pStableCells(1) = 1;
    goodCells = cell(1,nSessions);
    for s=2:nSessions
        load(fullfile(md(s).Location,'graphData_p.mat'),'A');
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
    keyboard;
    
end