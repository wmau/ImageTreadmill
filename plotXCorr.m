function plotXCorr(R,source,sink,lags,labels,controls)
%
%
%

%%
    %Set up error bars. 
    trialErrBar = zeros([size(R.trialShuffled{source,sink}.mu),2]);
    trialErrBar(1,:,1) = abs(R.trialShuffled{source,sink}.upper' - R.trialShuffled{source,sink}.mu');
    trialErrBar(1,:,2) = abs(R.trialShuffled{source,sink}.lower' - R.trialShuffled{source,sink}.mu'); 
%    errBar = smooth(R.shuffled{source,sink}.std)';
    trialLine.col = {'b'};
    
    timeErrBar = zeros([size(R.timeShuffled{source,sink}.mu),2]);
    timeErrBar(1,:,1) = abs(R.timeShuffled{source,sink}.upper' - R.timeShuffled{source,sink}.mu');
    timeErrBar(1,:,2) = abs(R.timeShuffled{source,sink}.lower' - R.timeShuffled{source,sink}.mu'); 
    timeLine.col = {'c'};
    
%     RErrBar = zeros([size(R.timeShuffled{source,sink}.mu),2]);
%     RErrBar(1,:,1) = abs(R.CI{source,sink}(1,:)' - R.smoothed{source,sink}');
%     RErrBar(1,:,2) = abs(R.CI{source,sink}(2,:)' - R.smoothed{source,sink}'); 
    Rline.col = {'k'};
    
    %Plot.
    %figure;
    hold on;
    
    if controls
        mseb(lags,R.trialShuffled{source,sink}.mu',trialErrBar,trialLine,1);    %Shuffled.
        mseb(lags,R.timeShuffled{source,sink}.mu',timeErrBar,timeLine,1);
    end
        
    %plot(lags,R.smoothed{source,sink},'k','linewidth',3);     %CCG.
    mseb(lags,R.smoothed{source,sink},R.CI{source,sink},Rline,1);
    sigplot = R.smoothed{source,sink}.*single(R.sig{source,sink});
    sigplot(sigplot==0) = nan;
    plot(lags,sigplot,'g','linewidth',3);                      %Significant.
    
    %Peak.  
    peak = round(nanmean(lags(R.sig{source,sink})),3);
    
    if isnan(peak)
        [~,peakind] = max(R.smoothed{source,sink});
        peak = lags(peakind); 
    end
    
    %Draw line showing average time at which xcorr is significant. 
    YLim = get(gca,'ylim');
    l = line([peak peak],YLim,'color','k','linestyle','--','linewidth',2);
    label(l,[num2str(peak*1000),' ms'],'location','top');
    
    if labels
        xlabel('Lags [s]');
        ylabel('Trial Averaged Cross-Correlation');
        title(['Neuron ',num2str(source),' x Neuron ',num2str(sink)]);
    end
    ylim(YLim);
end