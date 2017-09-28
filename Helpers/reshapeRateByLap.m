function X = reshapeRateByLap(md,varargin)
%X = reshapeRateByLap(md,varargin)
%
%   Reshapes ratebylap matrix into an observation x cell matrix as a
%   predictor matrix for Bayesian modeling. 
%

%% Parse inputs. 
    cd(md.Location);
    load('TimeCells.mat','ratebylap','TodayTreadmillLog');
    complete = logical(TodayTreadmillLog.complete);         %Completed runs.
    nRuns = round(sum(complete)*.5);                        %Number of runs to use in model.
    randRuns = randsample(sum(complete),nRuns)';            %Subset of runs. 
    
    p = inputParser; 
    p.addRequired('md',@(x) isstruct(x)); 
    p.addParameter('neurons',getTimeCells(md),@(x) isnumeric(x));
    p.addParameter('runs',randRuns,@(x) isnumeric(x)); 

    p.parse(md,varargin{:});
    neurons = p.Results.neurons;
    runs = p.Results.runs; 
    
    %Flip vector into column vector. 
    if size(neurons,1) > size(neurons,2)
        neurons = neurons';
    end
%% Reshape ratebylap. 
    ratebylap = ratebylap(runs,:,neurons);                  %Only take specified runs.
    ratebylap(ratebylap > 0) = 1; 
    [nRuns,nBins,nNeurons] = size(ratebylap);              	%Size variables. 
    
    %Preallocate then reshape. 
    X = nan(nRuns*nBins,nNeurons);
    for n=1:nNeurons
        temp = ratebylap(:,:,n)';           %Get that neuron's raster. 
        X(:,n) = temp(:);                   %Flatten. 
    end
end