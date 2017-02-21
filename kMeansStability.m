function [accuracy,sAccuracy,p] = kMeansStability(mds,stabilityCriterion,predictor,distType)
%
%
%

%%
    STATS = PartitionStats(mds,stabilityCriterion,predictor);
    B = 1000;
    
    nAnimals = length(STATS.stable);
    sStats = []; usStats = [];
    for a=1:nAnimals
        nSample = min([length(STATS.stable{a}), length(STATS.unstable{a})]);
        
        sStats = [sStats; randsample(STATS.stable{a},nSample)];
        usStats = [usStats; randsample(STATS.unstable{a},nSample)];  
    end
        
    nStats = size(sStats,1);
    
    X = [sStats; usStats];
    
    bad = find(X==0);
    X(bad) = [];
         
    clusters = kmeans(X,2,'distance',distType);
    
    if mean(X(clusters==2)) < mean(X(clusters==1))
        s=1; u=2;
    else
        s=2; u=1;
    end
    
    %Define labels.
    stable = [  s*ones(nStats,1);...
                u*ones(nStats,1)];
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