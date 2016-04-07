function scatterTimeSpatialInfo(MD)
%
%
%

%%
    dateTitle = MD.Date;
    dateTitle(3:3:6) = '-';
    load(fullfile(MD.Location,'PlaceMaps.mat'),'SpatialI');
    load(fullfile(MD.Location,'TimeCells.mat'),'TimeCells');
    load(fullfile(MD.Location,'TemporalInfo.mat'),'I'); 
    
    [~,p] = corr(SpatialI(TimeCells)',I(TimeCells),'type','spearman'); 
    scatter(SpatialI(TimeCells)',I(TimeCells),[],'.'); lsline;
    xlabel('Spatial Info [bits/s]'); 
    ylabel('Temporal Info [bits/s]'); 
    title({dateTitle,['p=',num2str(p)]}); 
end