disp('Partitioning temporal information based on temporal stability.');
[~,sTimeTIaccuracy,sTimeTIshuffle,sTimeTIp] = ClassifyStability(fulldataset,'time','TI');

disp('Partitioning spatial information based on temporal stability.');
[~,sTimeSIaccuracy,sTimeSIshuffle,sTimeSIp] = ClassifyStability(fulldataset,'time','SI');

disp('Partitioning temporal information based on spatial stability.');
[~,sPlaceTIaccuracy,sPlaceTIshuffle,sPlaceTIp] = ClassifyStability(fulldataset,'place','TI');

disp('Partitioning spatial information based on spatial stability.');
[~,sPlaceSIaccuracy,sPlaceSIshuffle,sPlaceSIp] = ClassifyStability(fulldataset,'place','SI');

B = length(sTimeTIshuffle);
edges = [.2:.02:.8];

figure;
    subplot(2,2,1);
    histogram(sTimeTIshuffle,edges,'normalization','probability','edgecolor','none'); 
    xlim([.2 .8]);
    yLim = get(gca,'ylim');
    line([.5 .5],[0 yLim(2)],'color','k','linewidth',2,'linestyle','--');
    line([sTimeTIaccuracy sTimeTIaccuracy],[0 yLim(2)],'color','r','linewidth',2);
    ylabel(['p = ',num2str(sTimeTIp)]);
    title('Temporal Stability Categorization from TIs');
    
    subplot(2,2,2);
    histogram(sTimeSIshuffle,edges,'normalization','probability','edgecolor','none'); 
    xlim([.2 .8]);
    yLim = get(gca,'ylim');
    line([.5 .5],[0 yLim(2)],'color','k','linewidth',2,'linestyle','--');
    line([sTimeSIaccuracy sTimeSIaccuracy],[0 yLim(2)],'color','r','linewidth',2);
    ylabel(['p = ',num2str(sTimeSIp)]);
    title('Temporal Stability Categorization from SIs');
    
    subplot(2,2,3);
    histogram(sPlaceTIshuffle,edges,'normalization','probability','edgecolor','none'); 
    xlim([.2 .8]);
    yLim = get(gca,'ylim');
    line([.5 .5],[0 yLim(2)],'color','k','linewidth',2,'linestyle','--');
    line([sPlaceTIaccuracy sPlaceTIaccuracy],[0 yLim(2)],'color','r','linewidth',2);
    ylabel(['p = ',num2str(sPlaceTIp)]);
    title('Spatial Stability Categorization from TIs');
    
    subplot(2,2,4);
    histogram(sPlaceSIshuffle,edges,'normalization','probability','edgecolor','none'); 
    xlim([.2 .8]);
    yLim = get(gca,'ylim');
    line([.5 .5],[0 yLim(2)],'color','k','linewidth',2,'linestyle','--');
    line([sPlaceSIaccuracy sPlaceSIaccuracy],[0 yLim(2)],'color','r','linewidth',2);
    ylabel(['p = ',num2str(sPlaceSIp)]);
    title('Stable Stability Categorization from SIs');