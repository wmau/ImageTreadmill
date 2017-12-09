function [Mdl,accuracy,shuffle,p] = LogRegressStability(mds,stabilityCriterion,predictor)
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
    
%%                       
    Mdl = fitclinear(X,stable,'learner','logistic','ClassNames',[1 0],...
        'Cost',C,'kfold',5);
    accuracy = 1-kfoldLoss(Mdl);

%% 
    shuffle = zeros(B,1);
    p = ProgressBar(B);
    parfor i=1:B
        r = stable(randperm(length(stable)));
        rMdl = fitclinear(X,r,'learner','logistic','classnames',[1 0],...
            'cost',C,'kfold',5);
        shuffle(i) = 1-kfoldLoss(rMdl);
        
        p.progress;
    end
    p.stop;
    
    %p-value.
    p = sum(accuracy <= shuffle)/B;
end