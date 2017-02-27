function [Mdl,accuracy,shuffle,p] = ClassifyStability(mds,stabilityCriterion,predictor,krnl)
%[Mdl,accuracy,shuffle,p] = ClassifyStability(mds,stabilityCriterion,predictor)
%
%   Trains a support vector machine on classifying stability in time or
%   space based on temporal or spatial information. Next, it trains B SVMs
%   on the same problem, but with shuffled stability labels. Then it
%   performs 10-fold cross validation on the real and shuffled data. 
%
%   INPUTS
%       mds: Session entries.
%
%       stabilityCriterion: Criterion for whether a neuron is labeled
%       stable. Can be either 'time' or 'place.
%
%       predictor: Statistic to use for categorization. Can be either 'SI'
%       or 'TI'.
%
%   OUTPUTS
%       Mdl: Output model trained to labeled neurons as either stable or
%       unstable based on a predictor. 
%
%       accuracy: Accuracy based on 10-fold cross validation.
%
%       shuffle: Distribution of accuracies for models trained on shuffled
%       labels. 
%
%       p: P-value of accuracy relative to shuffle. 
%

%% Set up.
    STATS = PartitionStats(mds,stabilityCriterion,predictor);
    B = 500;
    
    %Concatenate statistics.
    sStats = cell2mat(STATS.stable');
    usStats = cell2mat(STATS.unstable');
    
    %FOR DOWNSAMPLING UNCOMMENT THESE.
%     sStats = randsample(sStats,165);
%     usStats = randsample(usStats,81);
    
    %Determine number of stable and unstable cells. 
    nStable = length(sStats); 
    nUnstable = length(usStats);
    
    %Get the lower number. This will be the sample size for each group. 
    %sampSize = round(min([nStable,nUnstable])*.7);
    
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
                
%% Train support vector machine.
    ho = .7;
    Mdl = fitcsvm(X,stable,'KernelFunction',krnl,'ClassNames',[1 0],...
        'Cost',C);
    CV = crossval(Mdl,'kfold',5);
    accuracy = 1-kfoldLoss(CV); 
    
%% Train support vector machine after shuffling stability labels. 
    shuffle = zeros(B,1);
    disp('Training models based on shuffled labels.');
    resolution = 2;
    updateInc = round(B/(100/resolution));
    p = ProgressBar(100/resolution);
    parfor i=1:B
        r = stable(randperm(length(stable)))
        rMdl = fitcsvm(X,r,'KernelFunction',krnl,'ClassNames',[1 0],...
            'Cost',C);
        rCV = crossval(rMdl,'kfold',5);
        shuffle(i) = 1-kfoldLoss(rCV); 
        
        if round(i/updateInc) == (i/updateInc), p.progress; end
    end
    p.stop;
    
    %p-value.
    p = sum(accuracy <= shuffle)/B;
    
end