function [TI,cellID] = msTempInfo(mds,neurons)
%
%
%

%%
    nSessions = length(mds);
    nNeurons = length(neurons);
    
    %Map cells. 
    mapMD = getMapMD(mds);
    map = msMatchCells(mds,neurons,true);
    cellID = map(:,1);
    
    TI = nan(nNeurons,nSessions);
    for s=1:nSessions
        cd(mds(s).Location); 
        load('TemporalInfo.mat','MI');
        
        for n=1:nNeurons
            if map(n,s)>0 
                TI(n,s) = MI(map(n,s)); 
            end
        end
    end
    
end