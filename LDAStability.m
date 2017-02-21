function [Mdl,accuracy,shuffle,p] = LDAStability(mds,stabilityCriterion,predictor)
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
    
    X = [   sStats;...
            usStats];
        
    %Define labels.
    stable = [  ones(nStable,1);...
                zeros(nUnstable,1)];
            
    [~,minority] = min([nStable nUnstable]);
    [~,majority] = max([nStable nUnstable]);
    C = zeros(2);
    if majority==2             %If there are more unstable...
        C(minority,majority) = 1;
        C(majority,minority) = nStable/nUnstable;
    elseif majority==1         %If there are more stable...
        C(minority,majority) = 1;
        C(majority,minority) = nUnstable/nStable;
    end
            
    Mdl = fitcdiscr(X,stable,'Cost',C);
    %CV = crossval(Mdl);
    accuracy = 1-resubLoss(Mdl);
    
    p = ProgressBar(B);
    parfor i=1:B
        r = stable(randperm(length(stable)))
        rMdl = fitcdiscr(X,r,'Cost',C);
        %rCV = crossval(rMdl); 
        %shuffle(i) = 1-kfoldLoss(rCV); 
        shuffle(i) = 1-resubLoss(rMdl);
        
        p.progress;
    end
    p.stop;
    
    %p-value.
    p = sum(accuracy <= shuffle)/B;
end