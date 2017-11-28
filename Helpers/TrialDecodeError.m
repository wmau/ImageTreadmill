function percentCorrect = TrialDecodeError(decodedTrial,realTrial,trialBlockLims)
%
%
%

%%
    realTrialBlocks = getTrialBlockNum(realTrial,trialBlockLims)'; 
    percentCorrect = sum((realTrialBlocks - decodedTrial)==0)/length(decodedTrial); 
end