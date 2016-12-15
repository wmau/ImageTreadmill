disp('Partitioning temporal information based on temporal stability.');
[~,sTimeTIaccuracy,sTimeTIshuffle,sTimeTIp] = ClassifyStability(fulldataset,'time','TI');

disp('Partitioning spatial information based on temporal stability.');
[~,sTimeSIaccuracy,sTimeSIshuffle,sTimeSIp] = ClassifyStability(fulldataset,'time','SI');

disp('Partitioning temporal information based on spatial stability.');
[~,sPlaceTIaccuracy,sPlaceTIshuffle,sPlaceTIp] = ClassifyStability(fulldataset,'place','TI');

disp('Partitioning spatial information based on spatial stability.');
[~,sPlaceSIaccuracy,sPlaceSIshuffle,sPlaceSIp] = ClassifyStability(fulldataset,'place','SI');

B = length(sTimeTIshuffle);
l = .05 * B; 
u = .95 * B;
sTimeTIshuffle = sort(sTimeTIshuffle);
sPlaceSIshuffle = sort(sPlaceSIshuffle);
sPlaceTIshuffle = sort(sPlaceTIshuffle);
sTimeSIshuffle = sort(sTimeSIshuffle); 

%%
    acc = [sTimeTIaccuracy sPlaceSIaccuracy sPlaceTIaccuracy sTimeSIaccuracy]';
    e = [   sTimeTIshuffle(l)-mean(sTimeTIshuffle), sTimeTIshuffle(u)-mean(sTimeTIshuffle);...
            sPlaceSIshuffle(l)-mean(sPlaceSIshuffle),sPlaceSIshuffle(u)-mean(sPlaceSIshuffle);...
            sPlaceTIshuffle(l)-mean(sPlaceTIshuffle),sPlaceTIshuffle(u)-mean(sPlaceTIshuffle);...
            sTimeSIshuffle(l)-mean(sTimeSIshuffle),sTimeSIshuffle(u)-mean(sTimeSIshuffle)];
    e = abs(e);
    x = 1:4; 
    y = [   mean(sTimeTIshuffle),...
            mean(sPlaceSIshuffle),...
            mean(sPlaceTIshuffle),...
            mean(sTimeSIshuffle)];
    figure; hold on;
    for i=1:4
        b(i) = bar(i,acc(i),.5);
    end
    [b([1,4]).FaceColor] = deal([0 .5 .5]);
    [b([2,3]).FaceColor] = deal([0.5765 0.4392 0.8588]);
    [b([2,4]).EdgeColor] = deal([0.5765 0.4392 0.8588]);
    [b([1,3]).EdgeColor] = deal([0 .5 .5]);
    [b.LineWidth] = deal(3);
    [h,p]=boundedline(x,y,e,'alpha');
    h.LineWidth = 2;
    h.Color = [.7 .7 .7];
    p.FaceColor = [.7 .7 .7];
    ylim([0.3 0.7]);
    set(gca,'tickdir','out','linewidth',2,'xtick',[]);
    ylabel('Accuracy');