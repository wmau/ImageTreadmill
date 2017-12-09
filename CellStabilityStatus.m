function [stable,exiting,entering] = CellStabilityStatus(md1,md2,cellType)
%[stable,exiting,entering] = CellStabilityStatus(md1,md2,cellType)
%
%   Determine which neurons are stable, entering, or exiting. 
%
%   INPUTS
%       md1 and md2: session entries.
%
%       cellType: 'time' or 'place' for time or place cells.
%
%   OUTPUTS
%       stable: indices of cells that are stable (tuning curve correlation
%       p < 0.01, Bonferroni corrected). 
%
%       exiting: indices of cells that are exiting (tuning curve
%       correlation p > 0.01, Bonferroni-corrected).
%
%       entering: indices of cells that were not coding in md1, but are in
%       md2. 
%

%% Get new time/place cells.
    switch cellType
        case 'time'
            entering = getNewTimeCells(md1,md2);
            cellGet = 'timecells';
        case 'place'
            entering = getNewPlaceCells(md1,md2);
            cellGet = 'placecells';
    end
    
%% Get stable/exiting cells. 
    neurons = AcquireTimePlaceCells(md1,cellGet);
    neurons = EliminateUncertainMatches([md1,md2],neurons); 
    
    %Do tuning curve correlation. 
    switch cellType
        case 'time'
            corrs = CorrTrdmllTrace(md1,md2,neurons);
        case 'place'
            corrs = CorrPlaceFields(md1,md2,neurons); 
    end
    
    %Stability criterion. 
    stblcrit = 0.01/length(neurons); 
    
    %Get stable and exiting cells. 
    stable = intersect(find(corrs(:,2) < stblcrit),neurons);
    exiting = intersect(find(corrs(:,2) > stblcrit | isnan(corrs(:,2))),neurons);
end