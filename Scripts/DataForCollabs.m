clear;
loadMD;
ALL = MD(292:309);
nSessions = length(ALL);

for i=1:nSessions
    cd(ALL(i).Location);
    load('TimeCells.mat','TodayTreadmillLog');
    load('TreadmillTraces.mat','DFDTTrdmll');
    load('Pos_align.mat','DFDTtrace','x_adj_cm','y_adj_cm'); 
    
    inds = TodayTreadmillLog.inds;
    complete = logical(TodayTreadmillLog.complete);
    
    onTM = inds(complete,:); 
    
    DATA(i).Animal = ALL(i).Animal;
    DATA(i).Date = ALL(i).Date; 
    DATA(i).TreadmillTraces = DFDTTrdmll;
    DATA(i).TimeCells = getTimeCells(ALL(i)); 
    DATA(i).OnTreadmill = onTM;
    DATA(i).Traces = DFDTtrace;
    DATA(i).x = x_adj_cm;
    DATA(i).y = y_adj_cm;
end