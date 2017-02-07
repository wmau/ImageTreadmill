function neurons = AcquireTimePlaceCells(md,cellGet)
%
%
%

%%
    PCcrit = .01;
%%
    switch cellGet
        case 'timecells',neurons = getTimeCells(md);
        case 'placecells',neurons = getPlaceCells(md,PCcrit);
        case 'dual',neurons = intersect(getTimeCells(md),...
                getPlaceCells(md,PCcrit));
    end
end