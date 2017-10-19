function TimeDecoderPosteriorProbabilityPlot(postProbs)
%
%
%

%% Prep. 
    %Make time vector. 
    nBins = size(postProbs,1); 
    t = linspace(0,10,nBins); 

    %Sort posterior probabilities for scaling color bar. 
    s = sort(postProbs(:)); 
    
    %Take the mean across runs. 
    PPMatrix = mean(postProbs,3); 
    
    %Get highest probability on each row. 
    [~,best] = max(PPMatrix); 

%% Plot. 
    figure; hold on;
        imagesc(t,t,PPMatrix);
            caxis([s(0.0125*length(s)) s(0.9875*length(s))]);
            colormap hot;
            axis equal; axis tight; 
            xlabel('Real time (s)');
            ylabel('Decoded time (s)');
            set(gca,'xtick',[0 5 10],'ytick',[0 5 10],'tickdir','out',...
                'ydir','reverse'); 
            c = colorbar; 
            set(c,'ytick',c.Limits);
        plot(t,t,'g','linestyle',':','linewidth',2);
        plot(t,t(best),'b','linewidth',2);
end