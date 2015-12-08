function [tuningcurve,shufflecurve,p,sigcurve] = TimeTuning(ratebylap,delays,T)
%TimeTuning(ratebylap,delay)
%
%   Takes a rate by lap matrix of size LxB (L=number of laps, B=number of
%   time bins) and determines whether there is significant time tuning
%   using a permutation method. A temporal tuning curve is calculated by
%   averaging calcium events across laps. Then that tuning curve is
%   circularly permuted 10,000 times and the tuning curve of the shuffled
%   data is calculated. Then we find the bins where the empirical data
%   exceeds that of the shuffled data 95% of the time. 
%
%   INPUTS
%       ratebylap: LxB matrix of calcium events in binned time lap by lap. 
%
%       delays: Output from getLapResponses, vector for each lap specifying
%       the length of delay.
%
%       T: Delay duratoin you want to examine time tuning for.
%
%   OUTPUTS
%       tuningcurve: Vector of size B, histogram of calcium events per time
%       bin. 
%
%       shufflecurve: IxB matrix (I=number of iterations) containing
%       shuffled tuning curves. Each row is a tuning curve, which is the
%       mean of permuted lap responses. 
%
%       p: Vector of size B, containing p-values for each bin equal to the
%       proportion of the shuffled data that has that bin exceeding that of
%       the empirical. 
%
%       sigcurve: Binary vector of size B, 1 where the region is
%       significantly different.
%   

%% Initialize
    ratebylap = ratebylap(delays==T,:); 
    B = 1000;                         %Number of iterations. 
    crit = 0.05;                       %Significance level.
    [nLaps,nBins] = size(ratebylap); 

    %Preallocate. 
    shufflelaps = nan(nLaps,nBins,B); 
    shufflecurve = nan(B,nBins); 
        
    %Tuning curve.
    tuningcurve = mean(ratebylap); 
    
%% Circular permutation test.
    %Construct matrix of temporal shifts. 
    shifts = randi([0,nBins],B,nLaps); 
    
    %For each permutation...
    for i=1:B
        
        %...for each lap...
        for thisLap=1:nLaps
            %...shift the response by a random amount. 
             shufflelaps(thisLap,:,i) = circshift(ratebylap(thisLap,:),[0,shifts(i,thisLap)]); 
        end
        
        %Then compute the mean of the shuffled time responses. 
        shufflecurve(i,:) = mean(shufflelaps(:,:,i));
       
    end
    
    %p-value is the fraction of times that the shuffled response turned out
    %to be higher than the empirical. 
    p = sum(shufflecurve > repmat(tuningcurve,B,1))/B; 
    
    %Parts of the tuning curve that are significantly different. 
    sigcurve = p < crit; 
end
