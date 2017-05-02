clear;
loadMD;
ALL = MD(292:309);
nSessions = length(ALL);

for i=1:nSessions
    cd(ALL(i).Location);
    load('TimeCells.mat','TodayTreadmillLog');
    load('TreadmillTraces.mat','ZTrdmll');
    
    inds = TodayTreadmillLog.inds;
    complete = logical(TodayTreadmillLog.complete);
    
    onTM = inds(complete,:); 
    
    DATA(i).Animal = ALL(i).Animal;
    DATA(i).Date = ALL(i).Date; 
    DATA(i).Traces = ZTrdmll;
    DATA(i).TimeCells = getTimeCells(ALL(i)); 
    DATA(i).OnTreadmill = onTM;
    
end