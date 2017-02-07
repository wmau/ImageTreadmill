function neurons = EliminateUncertainMatches(mds,neurons)
%
%
%

%%
    matches = msMatchCells(getMapMD(mds),mds,neurons,true);
    neurons = intersect(neurons,matches(:,1));
    
end