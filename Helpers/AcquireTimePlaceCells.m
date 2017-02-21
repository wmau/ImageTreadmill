function neurons = AcquireTimePlaceCells(md,cellGet)
%
%
%

%%
    PCcrit = .01;
%%
    switch cellGet
        case 'timecells'
            neurons = getTimeCells(md);
%             PCs = getPlaceCells(md,.01);
%             neurons = setdiff(neurons,PCs);
        case 'placecells'
            neurons = getPlaceCells(md,PCcrit);
%             TCs = getTimeCells(md); 
%             neurons = setdiff(neurons,TCs);
        case 'dual',neurons = intersect(getTimeCells(md),...
                getPlaceCells(md,PCcrit));
    end
end