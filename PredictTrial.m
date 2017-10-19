function [decodedTrial,postProbs] = PredictTrial(Mdl,X,varargin)
%
%
%

%%
    p = inputParser;
    p.addRequired('Mdl');
    p.addRequired('X'); 
    p.addParameter('nBins',40,@(x) isnumeric(x)); 
    p.addParameter('plotit',true,@(x) islogical(x)); 

    p.parse(Mdl,X,varargin{:});
    nBins = p.Results.nBins; 
    plotit = p.Results.plotit;
    
    nTestLaps = size(X,1) / nBins; 
%%
    postProbs = nan(nBins,nBins,nTestLaps); 
    decodedTrial = nan(nBins,nTestLaps);
    
    %Reshape X so that it's observations x neurons. 
    for l=1:nTestLaps
        %Indexing nonsense. 
        chunk = (l-1)*nBins+1:(l-1)*nBins+nBins;
        
        %Predict. 
        [decodedTrial(:,l),postProbs(:,:,l)] = predict(Mdl,X(chunk,:));
    end
    
%% Plot. 
    if plotit 
        TimeDecoderPosteriorProbabilityPlot(postProbs); 
    end
end