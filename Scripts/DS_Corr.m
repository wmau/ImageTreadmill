clear;
loadMD;

fulldataset = MD(292:309); 
nSessions = length(fulldataset); 
neuralActivity = 'transient';

[tDS,pDS,TI,SI,sigT,sigP] = deal(cell(nSessions,1));
for s=1:nSessions
    disp(['Analyzing ',fulldataset(s).Animal,' ',fulldataset(s).Date]);
    neurons = AcquireTimePlaceCells(fulldataset(s),'dual');
    
    [tDS{s},sigT{s}] = DecodeTimeWithSpatialLookup(fulldataset(s),...
        'plotit',false,'neurons',neurons,'neuralActivity',neuralActivity);    
    [pDS{s},sigP{s}]= DecodePlaceWithTemporalLookup(fulldataset(s),...
        'plotit',false,'neurons',neurons,'neuralActivity',neuralActivity);
    
    cd(fulldataset(s).Location);
    load('TemporalInfo.mat','MI');
    TI{s} = MI;
    load('SpatialInfo.mat','MI');
    SI{s} = MI; 
end

T = cell2mat(tDS);
P = cell2mat(pDS); 
TI = cell2mat(TI);
SI = cell2mat(SI); 

[r,p] = corr(P,T,'rows','complete');
figure;
s = scatter(P,T,50,'k');
alpha(s,.5);
lsline;
title(['R=',num2str(r),', p=',num2str(p)],'fontsize',15);
xlabel('Spatial DS','fontsize',15);
ylabel('Temporal DS','fontsize',15); 
set(gca,'tickdir','out','fontsize',12,'linewidth',4);
% xlim([0 1]);
% ylim([0 1]);