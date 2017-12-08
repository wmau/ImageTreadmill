function trialBlocks = getTrialBlockNum(laps,blockLims)
%trialBlocks = getTrialBlockNum(laps,blockLims)
%
%   Throws trials into bins of trial blocks. 
%

%%
    nBlocks = length(blockLims)-1;
    trialBlocks = nan(1,length(laps));
    for blockNum = 1:nBlocks
        trialBlocks(laps > blockLims(blockNum) & ...
            laps <= blockLims(blockNum+1)) = blockNum;
    end
    
end