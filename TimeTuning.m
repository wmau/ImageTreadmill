function [tuningcurve,shufflecurve,p,sigcurve,ci] = TimeTuning(ratebylap,delaysetting,complete,T)
%[tuningcurve,shufflecurve,p,sigcurve,ci] = TimeTuning(ratebylap,delays,T)
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
%       TodayTreadmillLog: 
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
%       ci: 2xB matrix containing the upper and lower confidence bounds
%       based on tuning curves from doing the circular shift bootstrap. 
%

%% Initialize
    %Select complete laps only. 
    ratebylap = ratebylap(delaysetting==T & complete,:); 
    ratebylap = ratebylap(:,~isnan(ratebylap(1,:)));
    
    B = 1000;                           %Number of iterations. 
    crit = 0.01;                        %Significance level.
    [nCompleteLaps,nBins] = size(ratebylap); 

    %Preallocate. 
    shufflelaps = nan(nCompleteLaps,nBins,B); 
    shufflecurve = nan(B,nBins); 
        
    %Tuning curve.
    tuningcurve = mean(ratebylap); 
    
%% Circular permutation test.
    %For each permutation...
    for i=1:B
        %Shuffle each trial. 
        shufflelaps(:,:,i) = permuteTime(ratebylap);
         
        %Then compute the mean of the shuffled time responses. 
        shufflecurve(i,:) = mean(shufflelaps(:,:,i));
       
    end
    
    %Get confidence intervals. 
    shufflecurve = sort(shufflecurve);
    idxci = round([0.975;0.025].*B);
    ci = shufflecurve(idxci,:);
    
    %p-value is the fraction of times that the shuffled response turned out
    %to be higher than the empirical. 
    p = sum(shufflecurve > repmat(tuningcurve,B,1))./B; 
    
    %Parts of the tuning curve that are significantly different. 
    sigcurve = p < crit; 
end
