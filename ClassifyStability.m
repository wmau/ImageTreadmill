function [Mdl,accuracy,shuffle,p] = ClassifyStability(mds,stabilityCriterion,predictor)
%[Mdl,accuracy,shuffle,p] = ClassifyStability(mds,stabilityCriterion,predictor)
%
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
    sampSize = min([nStable,nUnstable]);
    
    X = [   sStats(randsample(nStable,sampSize));...
            usStats(randsample(nUnstable,sampSize))];
    
    %Define labels.
    stable = [  ones(sampSize,1);...
                zeros(sampSize,1)];
                
%% Train support vector machine.
    Mdl = fitcsvm(X,stable,'KernelFunction','rbf','Standardize',true);
    CV = crossval(Mdl);
    accuracy = 1-kfoldLoss(CV); 
    
%% Train support vector machine after shuffling stability labels. 
    shuffle = zeros(B,1);
    disp('Training models based on shuffled labels.');
    p = ProgressBar(B);
    parfor i=1:B
        rMdl = fitcsvm(X,stable(randperm(length(stable))),'KernelFunction',...
            'rbf','Standardize',true);
        rCV = crossval(rMdl); 
        shuffle(i) = 1-kfoldLoss(rCV); 
        
        p.progress;
    end
    p.stop;
    
    %p-value.
    p = sum(accuracy < shuffle)/B;
    
end