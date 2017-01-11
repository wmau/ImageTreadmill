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
    B = 1000;
    
    %Concatenate statistics.
    sStats = cell2mat(STATS.stable');
    usStats = cell2mat(STATS.unstable');
    
    %Determine number of stable and unstable cells. 
    nStable = length(sStats); 
    nUnstable = length(usStats);
    
    %Get the lower number. This will be the sample size for each group. 
    sampSize = round(min([nStable,nUnstable])*.7);
    
    X = [   sStats(randsample(nStable,sampSize));...
            usStats(randsample(nUnstable,sampSize))];
    
    %Define labels.
    stable = [  ones(sampSize,1);...
                zeros(sampSize,1)];
                
%% Train support vector machine.
    Mdl = fitcsvm(X,stable,'KernelFunction',krnl,'Standardize',true);
    CV = crossval(Mdl);
    accuracy = 1-kfoldLoss(CV); 
    
%% Train support vector machine after shuffling stability labels. 
    shuffle = zeros(B,1);
    disp('Training models based on shuffled labels.');
    p = ProgressBar(B);
    parfor i=1:B
        rMdl = fitcsvm(X,stable(randperm(length(stable))),'KernelFunction',...
            krnl,'Standardize',true);
        rCV = crossval(rMdl); 
        shuffle(i) = 1-kfoldLoss(rCV); 
        
        p.progress;
    end
    p.stop;
    
    %p-value.
    p = sum(accuracy <= shuffle)/B;
    
end