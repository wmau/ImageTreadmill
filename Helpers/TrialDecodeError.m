function decodeError = TrialDecodeError(decodedTrial,realTrial,trialBlockLims)
%
%
%

%%
    realTrialBlocks = getTrialBlockNum(realTrial,trialBlockLims)'; 
    decodeError = abs(realTrialBlocks - decodedTrial); 
end