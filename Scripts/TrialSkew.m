clear;
loadMD;

cellType = 'placecells';
rasterType = 'place';

switch cellType
    case 'timecells'
        c = [0 .5 .5];
    case 'placecells'
        c = [.58 .44 .86];
end

B = 100;
fulldataset = MD(292:309); 

nSessions = length(fulldataset); 
[shuffleExample,peaks,skewness] = deal(cell(nSessions,1));
for s=1:nSessions
    skewness{s} = getAllSkewnesses(fulldataset(s),'cellType',cellType,...
        'rasterType',rasterType); 
    %shuffleExample{s} = getAllSkewnesses(fulldataset(s),'shuffle',true);
    [~,peaks{s}] = getTimePeak(fulldataset(s));
    
    %cd(fulldataset(s).Location);
    %load('TemporalInfo.mat','MI');
    %ti{s} = MI;
end
SKEWS = cell2mat(skewness);
%SHUFFLE = cell2mat(shuffleExample);
%scatter(SKEWS,PSKEWS,'.');

figure; hold on;
h = histogram(SKEWS,'edgecolor','none','normalization',...
    'probability','binwidth',0.02,'facecolor',c);
set(gca,'tickdir','out','linewidth',4,'fontsize',12);
xlabel('Within-session trial bias','fontsize',15);

shuffle = cell(nSessions,1);
shuffleSD = nan(B,1);
for i=1:B
    for s=1:nSessions
        shuffle{s} = getAllSkewnesses(fulldataset(s),'shuffle',true,...
            'cellType',cellType,'rasterType',rasterType);
    end
    
    if i==1
        histogram(cell2mat(shuffle),'binwidth',0.02,'facecolor',[.7 .7 .7],...
            'edgecolor','none','normalization','probability');
        ylabel(['Proportion of ',rasterType,' cells'],'fontsize',15); 
        
        keyboard;
    end
    shuffleSD(i) = nanstd(cell2mat(shuffle));
end
SD = nanstd(PSKEWS);

histogram(shuffleSD,'edgecolor','none','normalization','probability')
hold on
line([SD SD],[0 0.16],'color','r')
set(gca,'tickdir','out');
xlabel('Standard Deviation','fontsize',15);
ylabel('Proportion','fontsize',15);