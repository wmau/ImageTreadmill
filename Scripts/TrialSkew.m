clear;
loadMD;

B = 1000;
fulldataset = MD(292:309); 

nSessions = length(fulldataset); 
[shuffleExample,peaks,skewness] = deal(cell(nSessions,1));
for s=1:nSessions
    skewness{s} = getAllSkewnesses(fulldataset(s)); 
    shuffleExample{s} = getAllSkewnesses(fulldataset(s),'shuffle',true);
    [~,peaks{s}] = getTimePeak(fulldataset(s));
    
    cd(fulldataset(s).Location);
    load('TemporalInfo.mat','MI');
    ti{s} = MI;
end
SKEWS = cell2mat(skewness);
SHUFFLE = cell2mat(shuffleExample);

histogram(SKEWS,'edgecolor','none','normalization','probability');
set(gca,'tickdir','out','linewidth',4,'fontsize',12);
xlabel('Trial skewness','fontsize',15);
ylabel('Proportion','fontsize',15);

shuffle = cell(nSessions,1);
shuffleSD = nan(B,1);
for i=1:B
    for s=1:nSessions
        shuffle{s} = getAllSkewnesses(fulldataset(s),'shuffle',true);
    end
    
    shuffleSD(i) = nanstd(cell2mat(shuffle));
end
SD = nanstd(SKEWS);

histogram(shuffleSD,'edgecolor','none','normalization','probability')
hold on
line([SD SD],[0 0.16],'color','r')
set(gca,'tickdir','out');
xlabel('Standard Deviation','fontsize',15);
ylabel('Proportion','fontsize',15);