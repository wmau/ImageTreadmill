clear;
loadMD;
fulldataset = MD(292:309);
krnl = 'gaussian';

saveBool = true;
folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures';
timeFileName = fullfile(folder,'SVM time');
placeFileName = fullfile(folder,'SVM place');

if saveBool
    c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');
    
    if ~strcmp(c,'y')
        saveBool = false;
    end
end

disp('Classifying temporal stability based on temporal information.');
[~,sTimeTIaccuracy,sTimeTIshuffle,sTimeTIp] = ClassifyStability(fulldataset,'time','ti',krnl);

disp('Classifying temporal stability based on spatial information.');
[~,sTimeSIaccuracy,sTimeSIshuffle,sTimeSIp] = ClassifyStability(fulldataset,'time','si',krnl);

disp('Classifying temporal stability based on transient frequency.');
[~,sTimeFRaccuracy,sTimeFRshuffle,sTimeFRp] = ClassifyStability(fulldataset,'time','FR',krnl);

disp('Classifying spatial stability based on spatial information.');
[~,sPlaceSIaccuracy,sPlaceSIshuffle,sPlaceSIp] = ClassifyStability(fulldataset,'place','si',krnl);

disp('Classifying spatial stability based on temporal information.');
[~,sPlaceTIaccuracy,sPlaceTIshuffle,sPlaceTIp] = ClassifyStability(fulldataset,'place','ti',krnl);

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
    shuffle = [sTimeTIshuffle sTimeSIshuffle sTimeFRshuffle];
    grps = [zeros(B,1) ones(B,1) 2*ones(B,1)];
%     e = [   sTimeTIshuffle(l)-mean(sTimeTIshuffle), sTimeTIshuffle(u)-mean(sTimeTIshuffle);...
%             sTimeSIshuffle(l)-mean(sTimeSIshuffle),sTimeSIshuffle(u)-mean(sTimeSIshuffle);...
%             sTimeFRshuffle(l)-mean(sTimeFRshuffle),sTimeFRshuffle(u)-mean(sTimeFRshuffle)];
%    e = abs(e);
%    x = 1:3; 
    figure('Position',[680 585 360 390]); hold on;
    for i=1:3
        s(i) = scatter(i,acc(i),100,'filled');
    end
    s(1).CData = timecolor;
    s(2).CData = spacecolor;
    s(3).CData = frcolor;
    b = boxplot(shuffle(:),grps(:),'color',timecolor,'symbol','d',...
        'labels',{'TI','SI','Activity rate'});
    ylim([0.4 0.75]);
    set(gca,'tickdir','out','fontsize',12,'linewidth',4)
    ylabel('Accuracy');
    title(['p = ',num2str(sTimeTIp),', ',num2str(sTimeSIp),', ',num2str(sTimeFRp)]);
    if saveBool
        print(timeFileName,'-dpdf');
    end
    
%% Plot place stability.
    acc = [sPlaceTIaccuracy sPlaceSIaccuracy sPlaceFRaccuracy]';
    shuffle = [sPlaceTIshuffle sPlaceSIshuffle sPlaceFRshuffle];
    grps = [zeros(B,1) ones(B,1) 2*ones(B,1)];

    figure('Position',[680 585 360 390]); hold on;
    for i=1:3
        s(i) = scatter(i,acc(i),100,'filled');
    end
    s(1).CData = timecolor;
    s(2).CData = spacecolor;
    s(3).CData = frcolor;
    b = boxplot(shuffle(:),grps(:),'color',spacecolor,'symbol','d',...
        'labels',{'TI','SI','Activity rate'});
    ylim([0.4 0.75]);
    set(gca,'tickdir','out','fontsize',12,'linewidth',4)
    ylabel('Accuracy');
    title(['p = ',num2str(sPlaceTIp),', ',num2str(sPlaceSIp),', ',num2str(sPlaceFRp)]);
    if saveBool
        print(timeFileName,'-dpdf');
    end