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
    randRuns = randsample(1:sum(complete),nRuns)';            %Subset of runs. 
    
    p = inputParser; 
    p.addRequired('md',@(x) isstruct(x)); 
    p.addParameter('neurons',getTimeCells(md),@(x) isnumeric(x));
    p.addParameter('runs',randRuns,@(x) isnumeric(x));
    p.addParameter('shuffle',false,@(x) islogical(x)); 
    p.addParameter('collapseTrials',false,@(x) islogical(x)); 

    p.parse(md,varargin{:});
    neurons = p.Results.neurons;
    runs = p.Results.runs; 
    shuffle = p.Results.shuffle;
    collapseTrials = p.Results.collapseTrials; 
    
    %Flip vector into column vector. 
    if size(neurons,1) > size(neurons,2)
        neurons = neurons';
    end
%% Reshape ratebylap. 
    ratebylap = ratebylap(complete,:,:);
    try
    ratebylap = ratebylap(runs,:,neurons);                  %Only take specified runs.
    catch, keyboard; end
    ratebylap(ratebylap > 0) = 1; 
    
    [nRuns,nBins,nNeurons] = size(ratebylap);              	%Size variables. 
    
    if shuffle 
        ratebylap = ratebylap(:,:,randperm(nNeurons));
    end
    
    %Preallocate then reshape. 
    if collapseTrials
        X = nan(nRuns,nNeurons);
        
        for n=1:nNeurons
            temp = ratebylap(:,:,n);
            X(:,n) = sum(temp,2);
        end
    else
        X = nan(nRuns*nBins,nNeurons);
        for n=1:nNeurons
            temp = ratebylap(:,:,n)';           %Get that neuron's raster. 
            X(:,n) = temp(:);                   %Flatten. 
        end
    end
end