loadMD;
fulldataset = MD(292:309);
krnl = 'linear';

disp('Classifying temporal stability based on temporal information.');
[~,sTimeTIaccuracy,sTimeTIshuffle,sTimeTIp] = ClassifyStability(fulldataset,'time','TI',krnl);

disp('Classifying temporal stability based on spatial information.');
[~,sTimeSIaccuracy,sTimeSIshuffle,sTimeSIp] = ClassifyStability(fulldataset,'time','SI',krnl);

disp('Classifying temporal stability based on transient frequency.');
[~,sTimeFRaccuracy,sTimeFRshuffle,sTimeFRp] = ClassifyStability(fulldataset,'time','FR',krnl);

disp('Classifying spatial stability based on spatial information.');
[~,sPlaceSIaccuracy,sPlaceSIshuffle,sPlaceSIp] = ClassifyStability(fulldataset,'place','SI',krnl);

disp('Classifying spatial stability based on temporal information.');
[~,sPlaceTIaccuracy,sPlaceTIshuffle,sPlaceTIp] = ClassifyStability(fulldataset,'place','TI',krnl);

disp('Classifying spatial stability based on transient frequency.');
[~,sPlaceFRaccuracy,sPlaceFRshuffle,sPlaceFRp] = ClassifyStability(fulldataset,'time','FR',krnl);

timecolor = [0 .5 .5];
spacecolor = [0.5765 0.4392 0.8588];
frcolor = [.7 .7 .7];

B = length(sTimeTIshuffle);
l = .05 * B; 
u = .95 * B;
sTimeTIshuffle = sort(sTimeTIshuffle);
sTimeSIshuffle = sort(sTimeSIshuffle); 
sTimeFRshuffle = sort(sTimeFRshuffle);
sPlaceSIshuffle = sort(sPlaceSIshuffle);
sPlaceTIshuffle = sort(sPlaceTIshuffle);
sPlaceFRshuffle = sort(sPlaceFRshuffle);

%% Plot time stability.
    acc = [sTimeTIaccuracy sTimeSIaccuracy sTimeFRaccuracy]';
    e = [   sTimeTIshuffle(l)-mean(sTimeTIshuffle), sTimeTIshuffle(u)-mean(sTimeTIshuffle);...
            sTimeSIshuffle(l)-mean(sTimeSIshuffle),sTimeSIshuffle(u)-mean(sTimeSIshuffle);...
            sTimeFRshuffle(l)-mean(sTimeFRshuffle),sTimeFRshuffle(u)-mean(sTimeFRshuffle)];
    e = abs(e);
    x = 1:3; 
    y = [   mean(sTimeTIshuffle),...
            mean(sTimeSIshuffle),...
            mean(sTimeFRshuffle)];
    figure; hold on;
    for i=1:3
        b(i) = bar(i,acc(i),.5);
    end
    [b(1:3).FaceColor] = deal(timecolor);
    b(2).EdgeColor = spacecolor;
    b(3).EdgeColor = frcolor;
    [b.LineWidth] = deal(3);
    [h,p]=boundedline(x,y,e,'alpha');
    h.LineWidth = 2;
    h.Color = 'b';
    p.FaceColor = 'b';
    ylim([0.35 0.6]);
    set(gca,'tickdir','out','linewidth',2);
    ylabel('Accuracy');
    xticklabels({'Temporal Stability | TI','Temporal Stability | SI','Temporal Stability | TR'})
    
%% Plot place stability.
    acc = [sPlaceSIaccuracy sPlaceTIaccuracy sPlaceFRaccuracy]';
    e = [   sPlaceSIshuffle(l)-mean(sPlaceSIshuffle), sPlaceSIshuffle(u)-mean(sPlaceSIshuffle);...
            sPlaceTIshuffle(l)-mean(sPlaceTIshuffle),sPlaceTIshuffle(u)-mean(sPlaceTIshuffle);...
            sPlaceFRshuffle(l)-mean(sPlaceFRshuffle),sPlaceFRshuffle(u)-mean(sPlaceFRshuffle)];
    e = abs(e);
    x = 1:3; 
    y = [   mean(sPlaceSIshuffle),...
            mean(sPlaceTIshuffle),...
            mean(sPlaceFRshuffle)];
    figure; hold on;
    for i=1:3
        b(i) = bar(i,acc(i),.5);
    end
    [b(1:3).FaceColor] = deal(spacecolor);
    b(2).EdgeColor = timecolor;
    b(3).EdgeColor = frcolor;
    [b.LineWidth] = deal(3);
    [h,p]=boundedline(x,y,e,'alpha');
    h.LineWidth = 2;
    h.Color = 'b';
    p.FaceColor = 'b';
    ylim([0.35 0.6]);
    set(gca,'tickdir','out','linewidth',2);
    ylabel('Accuracy');
    xticklabels({'Spatial Stability | SI','Spatial Stability | TI','Spatial Stability | TR'})
    xticklabels