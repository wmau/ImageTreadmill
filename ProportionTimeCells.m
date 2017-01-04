function [N,nTM,nTCs] = ProportionTimeCells(mds)
%
%
%

%%
    nSessions = length(mds); 
    [N,nTM,nTCs] = deal(zeros(1,nSessions));

    for s=1:nSessions
        cd(mds(s).Location);
        load('FinalOutput.mat','NumNeurons');
        load('TimeCells.mat','TimeCells');
        load('TemporalInfo.mat','sig');
        TimeCells = intersect(TimeCells,find(sig));
        
        N(s) = NumNeurons;
        nTM(s) = nNeuronsActiveonTM(mds(s)); 
        nTCs(s) = length(TimeCells);
    end
    
    A = round([mean(N),mean(nTM),mean(nTCs)]);
    I = round([mean(nTM),mean(nTCs)*ones(1,3)]);
     
    figure;
    venn(A,I,'facecolor',{'k','w','g'});
    axis equal; axis off; 
    
end