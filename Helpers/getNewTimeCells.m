function [newTCsS1,newTCsS2] = getNewTimeCells(base,comp)
%
%
%

%%
    DATA = CompileMultiSessionData([base,comp],{'timecells'});
    
    matches = msMatchMultiSessionCells([base,comp],DATA.timecells);
    
    %Get the cell identity in session one of all the time cells in
    %session 2. 
    newTCsS2 = matches(ismember(matches(:,2),DATA.timecells{2}),2);     %Cell numbers in session 2 of time cells in session 2.
    TCsfromS2 = matches(ismember(matches(:,2),DATA.timecells{2}),1);    %Cell numbers in session 1 of time cells in session 2. 
    TCsfromS1 = DATA.timecells{1};                                      %Cell numbers in session 1 of time cells in session 1. 
    
    [newTCsS1,ind] = setdiff(TCsfromS2,TCsfromS1);                      %Cell numbers in session 1 that became time cells in session 2. 
    newTCsS2 = newTCsS2(ind);                                           %Cell numbers in session 2 that became time cells in session 2. 
    
    %Erase neurons that didn't exist in session 1. 
    bad = newTCsS1 == 0 | isnan(newTCsS1);                              
    newTCsS1(bad) = [];
    newTCsS2(bad) = [];
end