dataset = MD(307:308);
nSessions = length(dataset); 
data = CompileMultiSessionData(dataset,{'timecells'});
map = msMatchMultiSessionCells(dataset,data.timecells);
map = map(all(map>0,2),:);
figure('Position',[680    60   420   915]);
for i=1:nSessions
    cd(dataset(i).Location); 
    load('TimeCells.mat','TodayTreadmillLog'); 
    nTrialBlocks = sum(TodayTreadmillLog.complete); 
    
    subplot(1,nSessions,i);
    if i==1
        skewness = getAllSkewnesses(dataset(i));
        [sorted,order] = sort(skewness(map(:,i)));
        trialDensityMap(dataset(i),'nTrialBlocks',nTrialBlocks,...
            'neurons',map(:,i),'order',order);
        notTC = find(double(isnan(sorted)),1,'first');
        line([1 nTrialBlocks],[notTC notTC],'color','r','linewidth',2);
    else
        trialDensityMap(dataset(i),'nTrialBlocks',nTrialBlocks,'order',order,...
            'neurons',map(:,i));
        line([1 nTrialBlocks],[notTC notTC],'color','r','linewidth',2);
    end
end