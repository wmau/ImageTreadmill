function [newPCsS1,newPCsS2] = getNewPlaceCells(base,comp)
%
%
%

%%
    DATA = CompileMultiSessionData([base,comp],{'placecells'});
    
    matches = msMatchMultiSessionCells([base,comp],DATA.placecells);
    
    %Get the cell identity in session one of all the time cells in
    %session 2. 
    newPCsS2 = matches(ismember(matches(:,2),DATA.placecells{2}),2);     %Cell numbers in session 2 of time cells in session 2.
    PCsfromS2 = matches(ismember(matches(:,2),DATA.placecells{2}),1);    %Cell numbers in session 1 of time cells in session 2. 
    PCsfromS1 = DATA.placecells{1};                                      %Cell numbers in session 1 of time cells in session 1. 
    
    [newPCsS1,ind] = setdiff(PCsfromS2,PCsfromS1);                      %Cell numbers in session 1 that became time cells in session 2. 
    newPCsS2 = newPCsS2(ind);                                           %Cell numbers in session 2 that became time cells in session 2. 
    
    %Erase neurons that didn't exist in session 1. 
    bad = newPCsS1 == 0 | isnan(newPCsS1);                              
    newPCsS1(bad) = [];
    newPCsS2(bad) = [];
end