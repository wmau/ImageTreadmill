function [tuningStatus,d] = PCRemap(ref,ssn)
%
%
%

%%
    ssns = [ref, ssn];
    nSessions = length(ssns);
    
    DATA = CompileMultiSessionData(ssns,{'placecells','placefieldcentroids'});
    PLACECELLS = DATA.placecells; 
    PFCentroids = DATA.placefieldcentroids; 
    nNeurons = length(PFCentroids); 
    
%% 
    mapMD = getMapMD(ref); 
    matchMat = msMatchCells(mapMD,ssns,PLACECELLS{1},true); 
    nPCs = size(matchMat,1); 
    
%% 
    [tuningStatus,d] = deal(nan(nNeurons,nSessions));
    for i=1:nPCs 
        n1 = matchMat(i,1); 
        
        for s=2:nSessions
            n2 = matchMat(i,s); 
            
            if isnan(n2) || n2==0, tuningStatus(n1,s) = nan; 
            %If not a place cell anymore, -1.
            elseif isempty(PFCentroids{s}{n2})
                tuningStatus(n1,s) = 0; 
            %If place cell, maybe centroid shifted. 
            else
                c1 = PFCentroids{1}{n1}; 
                c2 = PFCentroids{s}{n2};
                
                x1 = c1(1); y1 = c1(2); 
                x2 = c2(1); y2 = c2(2);
                
                d(n1,s-1) = sqrt((x2-x1)^2 + (y2-y1)^2); 
                
                if d(n1,s-1) > 9, tuningStatus(n1,s) = 0; 
                else, tuningStatus(n1,s) = 1; end
            end
            
        end
    end
end