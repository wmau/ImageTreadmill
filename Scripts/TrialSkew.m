clear;
loadMD;

B = 500;
fulldataset = MD(292:309); 

nSessions = length(fulldataset); 
[shuffleExample,peaks,Tskewness,Pskewness] = deal(cell(nSessions,1));
for s=1:nSessions
    Tskewness{s} = getAllSkewnesses(fulldataset(s),'cellType','timecells',...
        'rasterType','time'); 
    Pskewness{s} = getAllSkewnesses(fulldataset(s),'cellType','placecells',...
        'rasterType','place'); 
    %shuffleExample{s} = getAllSkewnesses(fulldataset(s),'shuffle',true);
    [~,peaks{s}] = getTimePeak(fulldataset(s));
    
    %cd(fulldataset(s).Location);
    %load('TemporalInfo.mat','MI');
    %ti{s} = MI;
end
TSKEWS = cell2mat(Tskewness);
PSKEWS = cell2mat(Pskewness);
%SHUFFLE = cell2mat(shuffleExample);
scatter(TSKEWS,PSKEWS,'.');

h = histogram(TSKEWS,'edgecolor','none','normalization',...
    'probability','binwidth',0.02);
set(gca,'tickdir','out','linewidth',4,'fontsize',12);
xlabel('Within-session trial bias','fontsize',15);
ylabel('Proportion','fontsize',15);

shuffle = cell(nSessions,1);
shuffleSD = nan(B,1);
for i=1:B
    for s=1:nSessions
        shuffle{s} = getAllSkewnesses(fulldataset(s),'shuffle',true,...
            'cellType','timecells','rasterType','time');
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