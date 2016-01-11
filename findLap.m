function lap = findLap(inds,treadmillEpochs)
%lap = findLap(inds,runEpochs)
%
%   Gets the lap number that correspond to the FT indices. Since we
%   sometimes specify "T", the delay duration of interest, we omit some
%   laps. Laps are also omitted if they are incomplete. The output of this
%   particular function states the lap number in the complete session
%   (i.e., without omitting laps). 
%
%   INPUTS
%       inds: Indices referencing FT for which we want the lap number for. 
%
%       runEpochs: Output from getTreadmillEpochs. Two column vector where
%       the first column is the index where the treadmill starts and the
%       second column is the index where the treadmill stops. 
%
%   OUTPUT
%       lap: Vector containing lap numbers for the indices input. NaN if
%       not during treadmill run. 
%

%% Find lap. 
    nInds = length(inds);
    lap = zeros(nInds,1); 
    
    for i=1:nInds
        try
            lap(i) = find(inds(i)>treadmillEpochs(:,1) & inds(i)<treadmillEpochs(:,2));
        catch
            lap(i) = nan;
        end          
    end
    
end