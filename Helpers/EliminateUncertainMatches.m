function neurons = EliminateUncertainMatches(mds,neurons)
%
%
%

%%
    matches = msMatchCells(mds,neurons,true);
    neurons = intersect(neurons,matches(:,1));
    
end