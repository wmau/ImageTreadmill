function TrialDecoderPosteriorProbabilityPlot(postProbs,trialBlockLims,testLaps)
%TrialDecoderPosteriorProbabilityPlot(postProbs,trialBlockLims,testLaps)
%
%   Plots the posterior probability matrix of a classifier output. Columsn
%   indicate data entries and rows are trial blocks. Color indicates
%   posterior probability. Green line is actual trial block identity. Blue
%   line is predicted trial block based on argmax posterior probability. 
%

%% Prep. 
    %Make trial vector. 
    nTestTrials = size(postProbs,1); 

    %Sort posterior probabilities for scaling color bar. 
    s = sort(postProbs(:)); 
    
    %Take the mean across runs. 
    PPMatrix = mean(postProbs,3); 
    
    %Get highest probability on each row. 
    [~,best] = max(PPMatrix,[],2); 
    
    %Get real trial block numbers. 
    realTrialBlocks = getTrialBlockNum(testLaps,trialBlockLims);
    trialBlocks = unique(realTrialBlocks); 

%% Plot. 
    figure; hold on;
        imagesc(1:nTestTrials,trialBlocks,PPMatrix');
            caxis([s(round(0.0125*length(s))) s(round(0.9875*length(s)))]);
            colormap hot;
            axis equal; axis tight; 
            xlabel('Test data entry');
            ylabel('Decoded trial block');
            set(gca,'tickdir','out','ydir','reverse'); 
            c = colorbar; 
            set(c,'ytick',c.Limits);
        plot(1:nTestTrials,realTrialBlocks,'g','linestyle',':','linewidth',2);
        plot(1:nTestTrials,trialBlocks(best),'b','linewidth',2);
end