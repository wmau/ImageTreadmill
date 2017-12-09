function [decodedTrial,postProbs,trialBlocks] = PredictTrial(Mdl,X,...
    trialBlockLims,testLaps,varargin)
%
%
%

%%
    p = inputParser;
    p.addRequired('Mdl');
    p.addRequired('X'); 
    p.addParameter('plotit',true,@(x) islogical(x)); 

    p.parse(Mdl,X,varargin{:});
    plotit = p.Results.plotit;
    
    trialBlocks = Mdl.ClassNames; 
    
    [decodedTrial,postProbs] = predict(Mdl,X); 
    
%% Plot. 
    if plotit 
        TrialDecoderPosteriorProbabilityPlot(postProbs,trialBlockLims,testLaps); 
    end
end