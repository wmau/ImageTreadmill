function [totalCells,TMcells,TCs] = ProportionTimeCells(mds)
%
%
%

%%
    nSessions = length(mds); 
    totalCells = zeros(1,nSessions);
    TMcells = zeros(1,nSessions);
    TCs = zeros(1,nSessions);
    
    for s=1:nSessions
        cd(mds(s).Location);
        load('Pos_align.mat','FT');
        load('TimeCells.mat','TimeCells');
        
        totalCells(s) = size(FT,1);
        TMcells(s) = nNeuronsActiveonTM(mds(s)); 
        TCs(s) = length(TimeCells);
    end
    
    figure;
    venn(round([mean(totalCells),mean(TMcells),mean(TCs)]),...
        round([mean(TMcells),mean(TCs)*ones(1,3)]),...
        'facecolor',{'k','w','g'});
    axis equal; axis off; 
    
end