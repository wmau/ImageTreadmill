function [accuracy,sAccuracy,p] = kMeansStability(mds,stabilityCriterion,predictor,distType)
%
%
%

%%
    STATS = PartitionStats(mds,stabilityCriterion,predictor);
    B = 1000;
    
    %Concatenate statistics.
    sStats = cell2mat(STATS.stable');
    usStats = cell2mat(STATS.unstable');
    
    %Determine number of stable and unstable cells. 
    nStable = length(sStats); 
    nUnstable = length(usStats);
    
    nSample = min([nStable, nUnstable]);

    X = [   randsample(sStats,nSample);...
            randsample(usStats,nSample)];
    
    bad = find(X==0);
    X(bad) = [];
         
    clusters = kmeans(X,2,'distance',distType);
    
    if median(X(clusters==2)) < median(X(clusters==1))
        s=1; u=2;
    else
        s=2; u=1;
    end
    
    %Define labels.
    stable = [  s*ones(nSample,1);...
                u*ones(nSample,1)];
    stable(bad) = [];
    
    accuracy = sum(stable==clusters)/length(stable);
    
    sAccuracy = nan(B,1);
    for i=1:B
        shuffleClusters = kmeans(X(randperm(length(stable))),2,'distance',distType);
        sAccuracy(i) = sum(stable==shuffleClusters)/length(stable);
    end
    
    p = sum(accuracy <= sAccuracy)/B;
    %keyboard;
end