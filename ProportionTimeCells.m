function [nTM,nTCs,nPCs] = ProportionTimeCells(mds)
%
%
%

%%
    nSessions = length(mds); 
    [nPCs,nTM,nTCs,PCandTM,PCandTC,all] = deal(zeros(1,nSessions));
    teal = [0 .5 .5];
    purple = [.58 .44 .86];

    PCcrit = .01;
    for s=1:nSessions
        cd(mds(s).Location);
                
        TimeCells = getTimeCells(mds(s));
        PlaceCells = getPlaceCells(mds(s),PCcrit);
        
        [nTM(s),active] = nNeuronsActiveonTM(mds(s)); 
        nTCs(s) = length(TimeCells);
        nPCs(s) = length(PlaceCells);
        
        PCandTM(s) = length(intersect(PlaceCells,active));
        PCandTC(s) = length(intersect(PlaceCells,TimeCells));
        all(s) = length(intersect(intersect(PlaceCells,TimeCells),active));
    end
    
    A = round([mean(nTM), mean(nTCs), mean(nPCs)]);
    I = round([mean(nTCs), mean(PCandTM), mean(PCandTC), mean(all)]);

    figure;
    venn(A,I,'facecolor',{'k',teal,purple});
    axis equal; axis off; 
    
end