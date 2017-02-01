function TCs = getTimeCells(md)
%
%
%

%%  
    cd(md.Location);
    load('TimeCells.mat','TimeCells');
    load('TemporalInfo.mat','sig'); 
    
    TCs = intersect(TimeCells,find(sig));
end