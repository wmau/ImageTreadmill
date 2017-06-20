clear;
loadMD;

fulldataset = MD(292:309); 
nSessions = length(fulldataset); 

[tDS,pDS] = deal(cell(nSessions,1));
for s=1:nSessions
    disp(['Analyzing ',fulldataset(s).Animal,' ',fulldataset(s).Date]);
    neurons = AcquireTimePlaceCells(fulldataset(s),'dual');
    
    tDS{s} = DecodeTimeWithSpatialLookup(fulldataset(s),'plotit',false,'neurons',neurons);    
    pDS{s} = DecodePlaceWithTemporalLookup(fulldataset(s),'plotit',false,'neurons',neurons);
end

T = cell2mat(tDS);
P = cell2mat(pDS); 

[r,p] = corr(P,T,'rows','complete');
figure;
scatter(P,T,40,[0 .5 .5],'filld'); 
title(['R=',num2str(r),', p=',num2str(p)],'fontsize',15);
xlabel('Spatial DS','fontsize',15);
ylabel('Temporal DS','fontsize',15); 
set(gca,'tickdir','out','fontsize',12,'linewidth',4);
xlim([0 1]);
ylim([0 1]);