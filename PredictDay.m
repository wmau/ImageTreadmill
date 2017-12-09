function [decodedDay,postProbs] = PredictDay(Mdl,X,varargin)
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
    
%% 
    [decodedDay,postProbs] = predict(Mdl,X); 
    
end