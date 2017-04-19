function TCs = getTimeCells(md)
%TCs = getTimeCells(md)
%
%   Gets time cells for that session. 
%

%% Main function.
    cd(md.Location);
    load('TimeCells.mat','TimeCells');
    load('TemporalInfo.mat','sig'); 
    
    TCs = intersect(TimeCells,find(sig));
    %TCs = TimeCells;
end