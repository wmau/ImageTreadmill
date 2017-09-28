function [decodedTime,postProbs] = PredictTime(Mdl,X,varargin)
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
    decodedTime = nan(nBins,nTestLaps);
    
    for l=1:nTestLaps
        chunk = (l-1)*nBins+1:(l-1)*nBins+nBins;
        
        [decodedTime(:,l),postProbs(:,:,l)] = predict(Mdl,X(chunk,:));
        postProbs(:,:,l) = postProbs(:,:,l)';
    end
    
    if plotit 
        figure;
        imagesc(mean(postProbs,3));
        axis equal;
    end
end