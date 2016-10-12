function plotXCorr(R,src,snk,lags,labels,controls)
%plotXCorr(R,src,snk,lags,labels,controls)
%
%

%% Check if src,snk is empty. 
    if isempty(R{src,snk})
        R = flipR(R,src,snk);
    end
    
%% Set up lines. 
    %Set up error bars. 
    trlErrBar = zeros([size(R{src,snk}.trlshff.mu),2]);
    trlErrBar(1,:,1) = abs(R{src,snk}.trlshff.upper' - R{src,snk}.trlshff.mu');
    trlErrBar(1,:,2) = abs(R{src,snk}.trlshff.lower' - R{src,snk}.trlshff.mu');     
    
    tErrBar = zeros([size(R{src,snk}.tshff.mu),2]);
    tErrBar(1,:,1) = abs(R{src,snk}.tshff.upper' - R{src,snk}.tshff.mu');
    tErrBar(1,:,2) = abs(R{src,snk}.tshff.lower' - R{src,snk}.tshff.mu'); 
    
    trlLine.col = {'b'};
    tLine.col = {'c'};
    Rline.col = {'k'};

%% Make plot. 
    figure;
    hold on;
    if controls
        mseb(lags,R{src,snk}.trlshff.mu',trlErrBar,trlLine,1);    %Shuffled.
        mseb(lags,R{src,snk}.tshff.mu',tErrBar,tLine,1);
    end
        
    mseb(lags,R{src,snk}.curve,R{src,snk}.CI,Rline,1);
    sigplot = R{src,snk}.curve.*single(R{src,snk}.sig);
    sigplot(sigplot==0) = nan;
    plot(lags,sigplot,'g','linewidth',3);                      %Significant.
    
    %Peak.  
    peak = round(nanmean(lags(R{src,snk}.sig)),3);
    
    if isnan(peak)
        [~,peakind] = max(R{src,snk}.curve);
        peak = lags(peakind); 
    end
    
    %Draw line showing average time at which xcorr is significant. 
    YLim = get(gca,'ylim');
    l = line([peak peak],YLim,'color','k','linestyle','--','linewidth',2);
    label(l,[num2str(peak*1000),' ms'],'location','top');
    
    if labels
        xlabel('Lags [s]');
        ylabel('Trial Averaged Cross-Correlation');
        title(['Neuron ',num2str(src),' x Neuron ',num2str(snk)]);
    end
    ylim(YLim);
end