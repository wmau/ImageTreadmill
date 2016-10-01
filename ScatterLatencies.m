function ScatterLatencies(md)
%
%
%

%%
    cd(md.Location);
    load('graphData_p.mat','A');
    load('TimeCells.mat','T','TodayTreadmillLog'); 
    
    complete = TodayTreadmillLog.complete;
    inds = TodayTreadmillLog.inds; 
    
    %Get treadmill run indices. 
    inds = inds(find(complete),:);  %Only completed runs. 
    inds(:,2) = inds(:,1) + 20*T-1; %Consistent length.   
    
    [~,targets] = find(A); 
    
    c=1; 
    nTargets = length(targets);
    cellSpread = cell(1,nTargets);
    treadmillSpread = cell(1,nTargets); 
    for n=targets'
        [~,cellSpread{c},treadmillSpread{c}] = SpreadRatio(md,A,n,'inds',inds);
        
        c = c+1;
    end
    
    figure; hold on;
    for n=1:nTargets
        for i=1:length(cellSpread{n})
            scatter(cellSpread{n}(i),treadmillSpread{n}(i),'ko','filled');
        end
    end
    
    XLim = get(gca,'xlim');
    YLim = get(gca,'ylim'); 
    Lim = max([XLim YLim]);
    plot([0:.1:Lim],[0:.1:Lim],'k');
    
    xlabel('Trigger-Target Latency MAD'); 
    ylabel('Treadmill-Target Latency MAD'); 
    set(gca,'ticklength',[0 0]);
    
end